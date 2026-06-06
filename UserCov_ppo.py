import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.distributions import Categorical
import math
import random
import time
import pandas as pd
import os

# 注意：假设 my_rtree 环境已配置
try:
    from my_rtree import RTree
except ImportError:
    class RTree:
        def get_height_3d(self, x, y): return 0.0
        def intersect3d(self, min_p, max_p): return False
        def load_from_file(self, path): pass

# --- 配置参数 ---
AREA_SIZE = 1000
B_TOTAL = 20e6
P_TX_DBM = 30.0
N0_DBM = -90.0
R_MIN_TH = 0.5e6
K_MAX = 50
FC_GHZ = 2.0
OMEGA = 0.85
N_STATIONS = 5
X_OFFSET = 835000
Y_OFFSET = 815000

class Rect3d:
    def __init__(self, minX=0, minY=0, minZ=0, maxX=0, maxY=0, maxZ=0):
        self.min = [minX, minY, minZ]
        self.max = [maxX, maxY, maxZ]

class CoverageEnv:
    def __init__(self, csv_path, rtree_path, num_users=600):
        self.num_users = num_users
        self.tree = RTree()
        if os.path.exists(rtree_path):
            self.tree.load_from_file(rtree_path)
        
        df_candidates = pd.read_csv(csv_path)
        # 确保使用列表下标作为 ID
        self.candidates = df_candidates.to_dict('records')
        self.num_candidates = len(self.candidates)
        self.users = self._generate_users()
        #保存到本地
        df=pd.DataFrame(self.users)
        df.to_csv('data/user_csv_ppo.csv',index=False,encoding='utf-8')
        print(f"用户坐标已保存到: {'data/user_csv_ppo.csv'}")

    def _generate_users(self):
        users = []
        for _ in range(int(self.num_users * 0.35)):
            r, theta = random.gauss(0, 50), random.uniform(0, 2 * math.pi)
            x = int(np.clip(250 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(250 + r * math.sin(theta), 0, AREA_SIZE))
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': True})
        
        for _ in range(int(self.num_users * 0.35)):
            r, theta = random.gauss(0, 80), random.uniform(0, 2 * math.pi)
            x = int(np.clip(750 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(750 + r * math.sin(theta), 0, AREA_SIZE))
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': True})
            
        while len(users) < self.num_users:
            x, y = random.uniform(0, AREA_SIZE), random.uniform(0, AREA_SIZE)
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': False})
        return users

    def get_path_loss(self, dist_3d, station_z, mode='LoS'):
        if mode == 'LoS':
            a, b = -28.0 - 20 * math.log10(FC_GHZ), 22.0
        else:
            a, b = 17.5 - 20 * math.log10(40 * math.pi * FC_GHZ / 3), 46 - 7 * math.log10(station_z)
        return -a + b * math.log10(max(dist_3d, 1e-3))

    def get_state(self, selected_indices):
        # 修正：返回 6 维向量
        if not selected_indices:
            return np.zeros(6, dtype=np.float32)
        
        sel_stations = [self.candidates[i] for i in selected_indices]
        mean_x = np.mean([s['x'] for s in sel_stations])
        mean_y = np.mean([s['y'] for s in sel_stations])
        
        state = [
            len(selected_indices) / N_STATIONS,
            (mean_x - 250) / AREA_SIZE,
            (mean_y - 250) / AREA_SIZE,
            (mean_x - 750) / AREA_SIZE,
            (mean_y - 750) / AREA_SIZE,
            len(selected_indices) * 0.1
        ]
        return np.array(state, dtype=np.float32)
    
    def evaluate(self, indices):
        selected_abs = [self.candidates[i] for i in indices]
        potential_links = [[] for _ in range(self.num_users)]
        
        for idx_in_indices, abs_row in enumerate(selected_abs):
            abs_pos = np.array([abs_row['x'], abs_row['y'], abs_row['z']])
            for u_idx, user in enumerate(self.users):
                u_pos = np.array([user['x'], user['y'], user['z']])
                dist_3d = np.linalg.norm(u_pos - abs_pos)
                
                mode = None
                if dist_3d <= abs_row['nlos_radius']:
                    rect = Rect3d(X_OFFSET + abs_pos[0], Y_OFFSET + abs_pos[1], abs_pos[2],
                                  X_OFFSET + u_pos[0], Y_OFFSET + u_pos[1], u_pos[2] + 1.5)
                    mode = 'NLoS' if self.tree.intersect3d(rect.min, rect.max) else 'LoS'
                elif dist_3d <= abs_row['los_radius']:
                    rect = Rect3d(X_OFFSET + abs_pos[0], Y_OFFSET + abs_pos[1], abs_pos[2],
                                  X_OFFSET + u_pos[0], Y_OFFSET + u_pos[1], u_pos[2] + 1.5)
                    if not self.tree.intersect3d(rect.min, rect.max):
                        mode = 'LoS'
                
                if mode:
                    pl = self.get_path_loss(dist_3d, abs_row['z'], mode)
                    snr = P_TX_DBM - pl - N0_DBM
                    # 使用列表中的位置索引作为临时标识
                    potential_links[u_idx].append({'id': indices[idx_in_indices], 'snr': snr})

        abs_loads = {i: 0 for i in indices}
        covered_count = 0
        hotspot_covered = 0
        total_hotspot = sum(1 for u in self.users if u['is_hotspot'])
        rates = []

        for u_idx, links in enumerate(potential_links):
            sorted_links = sorted(links, key=lambda x: x['snr'], reverse=True)
            for link in sorted_links:
                if abs_loads[link['id']] < K_MAX:
                    abs_loads[link['id']] += 1
                    bw_user = B_TOTAL / (N_STATIONS * abs_loads[link['id']])
                    snr_lin = 10**(link['snr']/10)
                    rate = bw_user * math.log2(1 + snr_lin)
                    if rate >= R_MIN_TH:
                        covered_count += 1
                        rates.append(rate)
                        if self.users[u_idx]['is_hotspot']: hotspot_covered += 1
                    break

        hs_ratio = hotspot_covered / total_hotspot if total_hotspot > 0 else 0
        total_ratio = covered_count / self.num_users
        affinity = OMEGA * hs_ratio + (1 - OMEGA) * total_ratio
        
        return {
            'affinity': affinity, 'hs_ratio': hs_ratio, 'total_ratio': total_ratio,
            'avg_rate': np.mean(rates) if rates else 0, 'min_rate': np.min(rates) if rates else 0,
            'avg_load': np.mean(list(abs_loads.values())) / K_MAX
        }

class ActorCritic(nn.Module):
    def __init__(self, n_inputs, n_outputs):
        super(ActorCritic, self).__init__()
        self.actor = nn.Sequential(
            nn.Linear(n_inputs, 128), nn.ReLU(),
            nn.Linear(128, 256), nn.ReLU(),
            nn.Linear(256, n_outputs), nn.Softmax(dim=-1)
        )
        self.critic = nn.Sequential(
            nn.Linear(n_inputs, 128), nn.ReLU(),
            nn.Linear(128, 1)
        )
    def forward(self, x):
        return self.actor(x), self.critic(x)

def train_ppo():
    csv_path = './data/coverage_results.csv'
    rtree_path = '/home/zoe/Documents/code/Q-view/mergeindex/indexes/ABS_selection.3idx'
    
    env = CoverageEnv(csv_path, rtree_path)
    # 修正：输入维度为 6
    input_dim = 6 
    output_dim = env.num_candidates
    model = ActorCritic(input_dim, output_dim)
    optimizer = optim.Adam(model.parameters(), lr=0.0003)
    
    max_episodes = 1000
    best_affinity = 0
    best_indices = []

    print(f"[*] 训练开始: 输入维度={input_dim}, 动作空间={output_dim}")
    t1 = time.time()
    
    for episode in range(max_episodes):
        selected_indices = []
        log_probs, values, rewards = [], [], []
        last_affinity = 0
        
        for step in range(N_STATIONS):
            state_np = env.get_state(selected_indices)
            state = torch.from_numpy(state_np).float()
            probs, val = model(state)
            
            mask = torch.ones(env.num_candidates)
            mask[selected_indices] = 0
            probs = probs * mask
            probs = probs / (probs.sum() + 1e-9)
            
            dist = Categorical(probs)
            action = dist.sample()
            
            selected_indices.append(action.item())
            log_probs.append(dist.log_prob(action))
            values.append(val)
            
            # Step Reward: 相比上一步的增量
            current_eval = env.evaluate(selected_indices)
            step_reward = current_eval['affinity'] - last_affinity
            rewards.append(step_reward)
            last_affinity = current_eval['affinity']

        # 更新网络
        returns = []
        gt = 0
        for r in reversed(rewards):
            gt = r + 0.95 * gt
            returns.insert(0, gt)
        
        returns = torch.tensor(returns, dtype=torch.float32)
        values = torch.cat(values).squeeze()
        log_probs = torch.stack(log_probs)
        
        advantages = returns - values.detach()
        actor_loss = -(log_probs * advantages).mean()
        critic_loss = nn.MSELoss()(values, returns)
        
        optimizer.zero_grad()
        (actor_loss + 0.5 * critic_loss).backward()
        optimizer.step()

        if last_affinity > best_affinity:
            best_affinity = last_affinity
            best_indices = selected_indices.copy()

        if episode % 20 == 0:
            print(f"Ep {episode:4d} | Best: {best_affinity:.4f} | Current: {last_affinity:.4f}")
            
    t2 = time.time()
    # 最终评估
    final_res = env.evaluate(best_indices)
    print("\n" + "="*60)
    print("      PPO REINFORCEMENT LEARNING FINAL RESULTS")
    print("="*60)
    print(f"最佳基站配置 (ID):      {best_indices}")
    print(f"热点区域覆盖率:         {final_res['hs_ratio'] * 100:.2f}%")
    print(f"总覆盖率:              {final_res['total_ratio']*100:.2f}%")
    print(f"热点区域平均速率:       {final_res['avg_rate'] / 1e6:.3f} Mbps")
    print(f"热点区域最低速率:       {final_res['min_rate'] / 1e6:.3f} Mbps")
    print(f"单机平均负载系数:       {final_res['avg_load'] * 100:.2f}%")
    print(f"运行时间 (1000 Ep):    {t2-t1:.2f}s")
    print("="*60)

if __name__ == "__main__":
    train_ppo()