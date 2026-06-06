import numpy as np
import pandas as pd
import random
import math
import copy
import time
import os

# --- 配置参数 ---
AREA_SIZE = 1000
B_TOTAL = 20e6        # 总带宽 20MHz
P_TX_DBM = 30.0       # 发射功率
N0_DBM = -90.0        # 噪声功率
R_MIN_TH = 0.5e6      # 速率阈值 0.5Mbps
K_MAX = 50            # 单基站最大负载用户数
FC_GHZ = 2.0          # 载波频率 2GHz
X_OFFSET = 835000     # 坐标偏移
Y_OFFSET = 815000
OMEGA = 0.85          # 亲和力权重 (热点覆盖权重)
N_STATIONS =  5      # 部署基站数量

# --- RTree 适配层 ---
try:
    from my_rtree import RTree
except ImportError:
    class RTree:
        def get_height_3d(self, x, y): return 0.0
        def intersect3d(self, min_p, max_p): return False
        def load_from_file(self, path): pass

class Rect3d:
    def __init__(self, minX, minY, minZ, maxX, maxY, maxZ):
        self.min = [minX, minY, minZ]
        self.max = [maxX, maxY, maxZ]

# --- 核心环境管理 ---
class CoverageEnv:
    def __init__(self, csv_path, rtree_path, num_users=600):
        self.num_users = num_users
        self.tree = RTree()
        if os.path.exists(rtree_path):
            self.tree.load_from_file(rtree_path)
        
        # 加载候选基站
        df = pd.read_csv(csv_path)
        # 确保包含 station_id，如果没有则用 index
        if 'station_id' not in df.columns:
            df['station_id'] = df.index
        self.candidates = df.to_dict('records')
        self.num_candidates = len(self.candidates)
        
        # 生成用户分布
        self.users = self._generate_users()
        #保存到本地
        df=pd.DataFrame(self.users)
        df.to_csv('data/user_csv_GA.csv',index=False,encoding='utf-8')
        print(f"用户坐标已保存到: {'data/user_csv_GA'}")

    def _generate_users(self):
        users = []
        # 热点 1: 30% 用户
        for _ in range(int(self.num_users * 0.35)):
            r, theta = random.gauss(0, 50), random.uniform(0, 2 * math.pi)
            x = int(np.clip(250 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(250 + r * math.sin(theta), 0, AREA_SIZE))
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': True})
        # 热点 2: 30% 用户
        for _ in range(int(self.num_users * 0.35)):
            r, theta = random.gauss(0, 80), random.uniform(0, 2 * math.pi)
            x = int(np.clip(750 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(750 + r * math.sin(theta), 0, AREA_SIZE))
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': True})
        # 散户
        while len(users) < self.num_users:
            x, y = int(random.uniform(0, AREA_SIZE)), int(random.uniform(0, AREA_SIZE))
            z = float(self.tree.get_height_3d(int(x) + X_OFFSET, int(y) + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': False})
        return users

    def get_path_loss(self, dist_3d, station_z, mode='LoS'):
        if mode == 'LoS':
            a, b = -28.0 - 20 * math.log10(FC_GHZ), 22.0
        else:
            a, b = 17.5 - 20 * math.log10(40 * math.pi * FC_GHZ / 3), 46 - 7 * math.log10(station_z)
        return -a + b * math.log10(max(dist_3d, 1e-3))

    def evaluate(self, indices):
        selected_bs = [self.candidates[i] for i in indices]
        potential_links = [[] for _ in range(self.num_users)]
        
        for bs in selected_bs:
            bs_pos = np.array([bs['x'], bs['y'], bs['z']])
            for u_idx, user in enumerate(self.users):
                u_pos = np.array([user['x'], user['y'], user['z']])
                dist_3d = np.linalg.norm(u_pos - bs_pos)
                
                mode = None
                if dist_3d <= bs['nlos_radius']:
                    rect = Rect3d(X_OFFSET+bs_pos[0], Y_OFFSET+bs_pos[1], bs_pos[2], X_OFFSET+u_pos[0], Y_OFFSET+u_pos[1], u_pos[2]+1.5)
                    mode = 'NLoS' if self.tree.intersect3d(rect.min, rect.max) else 'LoS'
                elif dist_3d <= bs['los_radius']:
                    rect = Rect3d(X_OFFSET+bs_pos[0], Y_OFFSET+bs_pos[1], bs_pos[2], X_OFFSET+u_pos[0], Y_OFFSET+u_pos[1], u_pos[2]+1.5)
                    if not self.tree.intersect3d(rect.min, rect.max): mode = 'LoS'
                
                if mode:
                    pl = self.get_path_loss(dist_3d, bs['z'], mode)
                    snr = P_TX_DBM - pl - N0_DBM
                    potential_links[u_idx].append({'id': bs['station_id'], 'snr': snr})

        bs_loads = {bs['station_id']: 0 for bs in selected_bs}
        covered_count, hotspot_covered = 0, 0
        total_hotspot = sum(1 for u in self.users if u['is_hotspot'])
        rates = []

        for u_idx, links in enumerate(potential_links):
            for link in sorted(links, key=lambda x: x['snr'], reverse=True):
                if bs_loads[link['id']] < K_MAX:
                    bs_loads[link['id']] += 1
                    bw_user = B_TOTAL / (len(indices) * bs_loads[link['id']])
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
            'indices': indices, 'affinity': affinity, 'hs_ratio': hs_ratio, 
            'total_ratio': total_ratio, 'avg_rate': np.mean(rates) if rates else 0,
            'min_rate': np.min(rates) if rates else 0, 'avg_load': np.mean(list(bs_loads.values())) / K_MAX
        }

# --- 遗传算法类 ---
class GeneticOptimizer:
    def __init__(self, env, pop_size=30, cx_pb=0.8, mut_pb=0.2):
        self.env = env
        self.pop_size = pop_size
        self.cx_pb = cx_pb
        self.mut_pb = mut_pb
        self.candidate_pool = list(range(env.num_candidates))

    def run(self, max_gen=200):
        # 1. 初始化种群
        population = [random.sample(self.candidate_pool, N_STATIONS) for _ in range(self.pop_size)]
        best_solution = None

        print(f"[*] Starting GA Optimization...")
        for gen in range(max_gen):
            # 2. 评估
            evaluated = [self.env.evaluate(ind) for ind in population]
            evaluated.sort(key=lambda x: x['affinity'], reverse=True)

            if not best_solution or evaluated[0]['affinity'] > best_solution['affinity']:
                best_solution = copy.deepcopy(evaluated[0])

            print(f"Gen {gen:3d} | Affinity: {best_solution['affinity']:.4f} | HS Cov: {best_solution['hs_ratio']*100:.2f}%")

            # 3. 选择与变异生成下一代
            next_pop = [best_solution['indices']] # 精英保留
            
            while len(next_pop) < self.pop_size:
                # 锦标赛选择
                p1 = self._tournament(evaluated)
                p2 = self._tournament(evaluated)
                
                # 交叉
                c1, c2 = self._cross(p1['indices'], p2['indices'])
                
                # 变异
                next_pop.append(self._mutate(c1))
                if len(next_pop) < self.pop_size: next_pop.append(self._mutate(c2))
            
            population = next_pop[:self.pop_size]

        return best_solution

    def _tournament(self, evaluated, k=3):
        return max(random.sample(evaluated, k), key=lambda x: x['affinity'])

    def _cross(self, p1, p2):
        if random.random() < self.cx_pb:
            point = random.randint(1, N_STATIONS - 1)
            # 保持唯一性交叉
            c1 = list(dict.fromkeys(p1[:point] + p2))[:N_STATIONS]
            c2 = list(dict.fromkeys(p2[:point] + p1))[:N_STATIONS]
            # 如果长度不足补齐
            while len(c1) < N_STATIONS: c1.append(random.choice(self.candidate_pool))
            while len(c2) < N_STATIONS: c2.append(random.choice(self.candidate_pool))
            return c1, c2
        return p1.copy(), p2.copy()

    def _mutate(self, ind):
        if random.random() < self.mut_pb:
            idx = random.randint(0, N_STATIONS - 1)
            ind[idx] = random.choice(self.candidate_pool)
        return ind

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    stations_csv = os.path.join(current_dir, 'data', 'coverage_results.csv')
    rtree_file = '/home/zoe/Documents/code/Q-view/mergeindex/indexes/ABS_selection.3idx'

    env = CoverageEnv(stations_csv, rtree_file)
    optimizer = GeneticOptimizer(env, pop_size=30)
    
    start_t = time.time()
    res = optimizer.run(max_gen=200)
    end_t = time.time()

    print("\n" + "="*55)
    print("      GENETIC ALGORITHM (GA) DEPLOYMENT RESULTS")
    print("="*55)
    print(f"部署 ID 集:            {res['indices']}")
    print(f"热点区域覆盖率:         {res['hs_ratio'] * 100:.2f}%")
    print(f"总覆盖率:              {res['total_ratio']*100:.2f}%")
    print(f"热点区域平均速率:       {res['avg_rate'] / 1e6:.3f} Mbps")
    print(f"热点区域最低速率:       {res['min_rate'] / 1e6:.3f} Mbps")
    print(f"单机平均负载系数:       {res['avg_load'] * 100:.2f}%")
    print(f"运行总时长:            {end_t - start_t:.2f}s")
    print("="*55)