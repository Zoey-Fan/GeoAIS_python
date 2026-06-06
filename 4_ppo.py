import time
import struct
import numpy as np
import pandas as pd
import os
import torch
import gymnasium as gym
from gymnasium import spaces
from stable_baselines3 import PPO
from scipy.spatial import KDTree

# 预计算比特计数表
_pop_range = np.arange(256, dtype=np.uint8)
POPCOUNT_TABLE = (_pop_range[:, None] & (1 << np.arange(8)) > 0).sum(1).astype(np.uint8)

class EnvironmentManager:
    """保持原有的数据加载逻辑不变"""
    def __init__(self, csv_path: str, bin_path: str):
        self.csv_path = csv_path
        self.bin_path = bin_path
        self.station_coords = None
        self.bitmaps = {}        
        self.individual_counts = {} 
        self.kd_tree = None
        self.all_station_ids = []
        self.BITMAP_SIZE = 125000  # 1000*1000 / 8

    def load_data(self):
        try:
            df = pd.read_csv(self.csv_path)
            self.station_coords = df.set_index('station_id')[['x', 'y', 'z']]
            self.all_station_ids = df['station_id'].tolist()
            self.kd_tree = KDTree(self.station_coords.values)
            self._load_binary_bitmaps()
        except Exception as e:
            print(f"[Error] {e}")

    def _load_binary_bitmaps(self):
        try:
            with open(self.bin_path, 'rb') as f:
                while True:
                    id_bytes = f.read(4)
                    if not id_bytes: break
                    station_id = struct.unpack('<i', id_bytes)[0]
                    bitmap_bytes = f.read(self.BITMAP_SIZE)
                    if len(bitmap_bytes) != self.BITMAP_SIZE: break
                    
                    if station_id in self.station_coords.index:
                        arr = np.frombuffer(bitmap_bytes, dtype=np.uint8).copy()
                        self.bitmaps[station_id] = arr         
                        self.individual_counts[station_id] = POPCOUNT_TABLE[arr].sum()
        except Exception as e:
            print(f"[Error] {e}")

# --- 强化学习环境定义 ---

class BaseStationEnv(gym.Env):
    """自定义强化学习环境"""
    def __init__(self, env_manager: EnvironmentManager, n_select: int):
        super(BaseStationEnv, self).__init__()
        self.mgr = env_manager
        self.n_select = n_select
        self.station_keys = list(self.mgr.bitmaps.keys())
        
        # 动作空间：选择基站的索引
        self.action_space = spaces.Discrete(len(self.station_keys))
        
        # 观测空间：为了计算效率，我们将1000x1000的位图下采样为32x32的特征图
        # 加上当前已选择的数量
        self.obs_dim = 32 * 32
        self.observation_space = spaces.Box(low=0, high=1, shape=(self.obs_dim + 1,), dtype=np.float32)
        
        self.reset()

    def _get_obs(self):
        # 将125000字节的位图转换为特征向量（约1000x1000像素下采样）
        # 这里使用简单步长采样代替复杂的池化操作以保证速度
        full_bits = np.unpackbits(self.current_union)  # 1000000维 (1000x1000)
        indices = np.linspace(0, len(full_bits)-1, self.obs_dim, dtype=int)
        downsampled = full_bits[indices].astype(np.float32)
        # 拼接步数特征（保证总维度1025）
        step_feature = np.array([self.steps / self.n_select], dtype=np.float32)
        return np.concatenate([downsampled, step_feature])

    def reset(self, seed=None, options=None):
        super().reset(seed=seed)
        self.current_union = np.zeros(self.mgr.BITMAP_SIZE, dtype=np.uint8)
        self.selected_indices = set()
        self.steps = 0
        self.total_coverage = 0
        return self._get_obs(), {}

    def step(self, action):
        action = int(action)
        station_id = self.station_keys[action]
        reward = 0
        terminated = False
        
        # 惩罚：重复选择同一个基站
        if action in self.selected_indices:
            reward = -20.0 
        else:
            self.selected_indices.add(action)
            candidate_bitmap = self.mgr.bitmaps[station_id]
            
            # 计算新增覆盖
            additional_bitmap = np.bitwise_and(candidate_bitmap, np.bitwise_not(self.current_union))
            additional_coverage = POPCOUNT_TABLE[additional_bitmap].sum()
            
            # 计算重叠惩罚 (新基站总面积 - 新增面积 = 重叠面积)
            overlap_area = self.mgr.individual_counts[station_id] - additional_coverage
            
            # 奖励函数设计：新增覆盖 - alpha * 重叠面积
            # 这里的 alpha 是权衡因子，可以根据具体需求调整
            reward = (0.7*additional_coverage / 1000) - (0.3 * overlap_area / 1000)
            
            # 更新状态
            self.current_union = np.bitwise_or(self.current_union, candidate_bitmap)
            self.total_coverage = POPCOUNT_TABLE[self.current_union].sum()
            
        self.steps += 1
        if self.steps >= self.n_select:
            terminated = True
            
        return self._get_obs(), reward, terminated, False, {}

# --- PPO 执行器 ---

class PPOSelector:
    def __init__(self, env_mgr: EnvironmentManager, n_stations: int):
        self.env_mgr = env_mgr
        self.n_stations = n_stations
        self.selected_ids = []
        self.total_coverage = 0
        self.overlap = 0.0
        self.run_time = 0.0
        
        # 创建环境
        self.gym_env = BaseStationEnv(env_mgr, n_stations)

    def calculate_additional_coverage(self, candidate_id, current_union):
        candidate_bitmap = self.env_mgr.bitmaps[candidate_id]
        additional_bitmap = np.bitwise_and(candidate_bitmap, np.bitwise_not(current_union))
        return POPCOUNT_TABLE[additional_bitmap].sum()

    def _rank_actions_by_policy(self, model, obs):
        obs_tensor = torch.as_tensor(obs, dtype=torch.float32, device=model.device).unsqueeze(0)
        with torch.no_grad():
            distribution = model.policy.get_distribution(obs_tensor)
            probs = distribution.distribution.probs.squeeze(0).cpu().numpy()
        return np.argsort(probs)[::-1]

    def _select_valid_action(self, model, obs, current_union):
        selected_indices = self.gym_env.selected_indices
        for action in self._rank_actions_by_policy(model, obs):
            action = int(action)
            if action not in selected_indices:
                return action

        # Fallback should only trigger if the policy ranking is unavailable or all
        # actions were already used. Use coverage gain to keep the final set valid.
        remaining_actions = [
            idx for idx in range(len(self.gym_env.station_keys))
            if idx not in selected_indices
        ]
        if not remaining_actions:
            return None

        return max(
            remaining_actions,
            key=lambda idx: self.calculate_additional_coverage(
                self.gym_env.station_keys[idx],
                current_union
            )
        )

    def run(self, train_steps=10000):
        """执行训练与推理"""
        print(f"\n[*] Starting PPO Reinforcement Learning for {self.n_stations} stations...")
        start_time = time.time()
        
        # 1. 训练模型 (对照实验建议至少训练一段时间，或者加载预训练模型)
        print(f"[-] Training PPO Agent for {train_steps} steps...")
        model = PPO(
            "MlpPolicy", 
            self.gym_env, 
            verbose=0, 
            learning_rate=0.0003,
            device="cpu"  # 强制CPU，避免GPU低效问题
        )
        model.learn(total_timesteps=train_steps)
        
        # 2. 推理执行
        print("[-] Executing optimized policy...")
        obs, _ = self.gym_env.reset()
        current_union = np.zeros(self.env_mgr.BITMAP_SIZE, dtype=np.uint8)
        
        for step in range(self.n_stations):
            action = self._select_valid_action(model, obs, current_union)
            if action is None:
                print(f"[Warning] No more stations to select (only {step} selected)")
                break

            station_id = self.gym_env.station_keys[action]
            
            # PPO's raw Discrete action may repeat a selected station. Inference
            # accepts the highest-probability action that is still feasible.
            self.selected_ids.append(station_id)
            
            # 计算新增覆盖（用于日志输出）
            candidate_bitmap = self.env_mgr.bitmaps[station_id]
            additional_bitmap = np.bitwise_and(candidate_bitmap, np.bitwise_not(current_union))
            max_additional = POPCOUNT_TABLE[additional_bitmap].sum()
            
            # 步进环境获取下一个观测
            obs, reward, terminated, _, _ = self.gym_env.step(action)
            
            # 更新最终统计
            current_union = np.bitwise_or(current_union, candidate_bitmap)
            self.total_coverage = POPCOUNT_TABLE[current_union].sum()
            
            sum_individual = sum(self.env_mgr.individual_counts[sid] for sid in self.selected_ids)
            self.overlap = (sum_individual - self.total_coverage) / self.total_coverage if self.total_coverage > 0 else 0.0
            
            print(f"Step {step+1}: Selected ID={station_id}, Reward={reward:.2f}, Total Coverage={self.total_coverage}, Overlap={self.overlap:.4f}")

        self.run_time = time.time() - start_time
        
        print("\n[=] PPO Optimization Finished.")
        print(f"Final Total Coverage: {self.total_coverage} pixels")
        print(f"Final Overlap Rate: {self.overlap:.4f}")
        print(f"Selected Station IDs: {self.selected_ids}")
        print(f"PPO Algorithm Run Time (Train+Run): {self.run_time:.6f} seconds")

if __name__ == "__main__":
    csv_path = 'data/coverage_results.csv'
    bin_path = 'data/coverage.bin'
    
    if not os.path.exists(csv_path) or not os.path.exists(bin_path):
        print("[!] Please provide data files.")
        exit()

    # 1. 加载数据
    env_mgr = EnvironmentManager(csv_path, bin_path)
    env_mgr.load_data()

    # 2. 配置参数
    n_stations = 50
    
    # 3. 运行 PPO 算法
    selector = PPOSelector(env_mgr, n_stations)
    # train_steps 决定了算法的“聪明”程度，对照实验可根据时间预算调整
    selector.run(train_steps=20000) 

    # 4. 输出汇总信息
    print("\n===== PPO Algorithm Summary =====")
    print(f"Number of Selected Stations: {n_stations}")
    print(f"Total Coverage: {selector.total_coverage} pixels")
    print(f"Overlap Rate: {selector.overlap:.4f}")
    print(f"Selected IDs: {selector.selected_ids}")
    print(f"Algorithm Run Time: {selector.run_time:.6f} s")
