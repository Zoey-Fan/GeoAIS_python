import numpy as np
import pandas as pd
import random
import math
import time
import os

# --- 注意：my_rtree 环境依赖 ---
try:
    from my_rtree import RTree
except ImportError:
    class RTree:
        def get_height_3d(self, x, y): return 0.0
        def intersect3d(self, min_p, max_p): return False
        def load_from_file(self, path): pass

# --- 配置参数 (与 AIS/RL/GA 代码严格一致) ---
AREA_SIZE = 1000
B_TOTAL = 20e6        # 总带宽 20MHz
P_TX_DBM = 30.0       # 发射功率
N0_DBM = -90.0        # 噪声功率
R_MIN_TH = 0.5e6      # 速率阈值 0.5Mbps
K_MAX = 50            # 单无人机最大服务用户量
FC_GHZ = 2.0          # 载波频率 2GHz
X_OFFSET = 835000
Y_OFFSET = 815000
OMEGA = 0.85          # 亲和力权重
N_STATIONS = 5        # 用户明确要求：生成5个基站

class Rect3d:
    def __init__(self, minX=0, minY=0, minZ=0, maxX=0, maxY=0, maxZ=0):
        self.min = [minX, minY, minZ]
        self.max = [maxX, maxY, maxZ]

# --- 环境管理 ---
class CoverageEnv:
    def __init__(self, csv_path, rtree_path, num_users=1000):
        self.num_users = num_users
        self.tree = RTree()
        if os.path.exists(rtree_path):
            self.tree.load_from_file(rtree_path)
        
        # 加载基站候选集
        if os.path.exists(csv_path):
            self.df_candidates = pd.read_csv(csv_path)
            self.candidates = self.df_candidates.to_dict('records')
            self.num_candidates = len(self.candidates)
        else:
            raise FileNotFoundError(f"Candidate file {csv_path} not found.")
            
        # 生成用户分布
        self.users = self._generate_users()
        #保存到本地
        df=pd.DataFrame(self.users)
        df.to_csv('data/user_csv_uniform.csv',index=False,encoding='utf-8')
        print(f"用户坐标已保存到: {'data/user_csv_uniform.csv'}")

    def _generate_users(self):
        users = []
        # 热点 1: 250, 250
        for _ in range(int(self.num_users * 0.35)):
            r = random.gauss(0, 50)
            theta = random.uniform(0, 2 * math.pi)
            x = int(np.clip(250 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(250 + r * math.sin(theta), 0, AREA_SIZE))
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': True})
        
        # 热点 2: 750, 750
        for _ in range(int(self.num_users * 0.35)):
            r = random.gauss(0, 80)
            theta = random.uniform(0, 2 * math.pi)
            x = int(np.clip(750 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(750 + r * math.sin(theta), 0, AREA_SIZE))
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': True})
            
        # 散户
        while len(users) < self.num_users:
            x = int(random.uniform(0, AREA_SIZE))
            y = int(random.uniform(0, AREA_SIZE))
            z = float(self.tree.get_height_3d(x + X_OFFSET, y + Y_OFFSET) or 0.0)
            users.append({'x': x, 'y': y, 'z': z, 'is_hotspot': False})
        return users

    def get_path_loss(self, dist_3d, station_z, mode='LoS'):
        if mode == 'LoS':
            a, b = -28.0 - 20 * math.log10(FC_GHZ), 22.0
        else:
            a, b = 17.5 - 20 * math.log10(40 * math.pi * FC_GHZ / 3), 46 - 7 * math.log10(station_z)
        return -a + b * math.log10(max(dist_3d, 1e-3))

    def evaluate(self, indices, n_stations=N_STATIONS):
        selected_abs = [self.candidates[i] for i in indices]
        potential_links = [[] for _ in range(self.num_users)]
        
        for abs_row in selected_abs:
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
                    potential_links[u_idx].append({'id': abs_row['station_id'], 'snr': snr})

        abs_loads = {self.candidates[i]['station_id']: 0 for i in indices}
        covered_count = 0
        hotspot_covered = 0
        total_hotspot = sum(1 for u in self.users if u['is_hotspot'])
        rates = []

        for u_idx, links in enumerate(potential_links):
            sorted_links = sorted(links, key=lambda x: x['snr'], reverse=True)
            for link in sorted_links:
                if abs_loads[link['id']] < K_MAX:
                    abs_loads[link['id']] += 1
                    bw_user = B_TOTAL / (n_stations * abs_loads[link['id']])
                    snr_lin = 10**(link['snr']/10)
                    rate = bw_user * math.log2(1 + snr_lin)
                    if rate >= R_MIN_TH:
                        covered_count += 1
                        rates.append(rate)
                        if self.users[u_idx]['is_hotspot']: hotspot_covered += 1
                    break

        hs_ratio = hotspot_covered / total_hotspot if total_hotspot > 0 else 0
        total_ratio = covered_count / self.num_users
        
        return {
            'indices': indices,
            'hs_ratio': hs_ratio, 
            'total_ratio': total_ratio,
            'avg_rate': np.mean(rates) if rates else 0, 
            'min_rate': np.min(rates) if rates else 0,
            'avg_load': np.mean(list(abs_loads.values())) / K_MAX
        }

def get_uniform_targets(n):
    """根据基站数量 n，生成在 1000x1000 区域内的均匀目标坐标"""
    if n == 5:
        return [(500, 500), (250, 250), (250, 750), (750, 250), (750, 750)]
    
    # 通用逻辑：尝试构建最接近平方根的网格
    rows = int(math.sqrt(n))
    cols = math.ceil(n / rows)
    
    targets = []
    for r in range(rows):
        for c in range(cols):
            if len(targets) < n:
                # 在每个网格单元的中心布放
                x = (c + 0.5) * (AREA_SIZE / cols)
                y = (r + 0.5) * (AREA_SIZE / rows)
                targets.append((x, y))
    return targets

def run_uniform_placement(env, n):
    """在范围内寻找 n 个均匀分布的候选基站"""
    targets = get_uniform_targets(n)
    selected_indices = []
    
    for tx, ty in targets:
        best_idx = -1
        min_dist = float('inf')
        for i, cand in enumerate(env.candidates):
            if i in selected_indices:
                continue
            dist = math.sqrt((cand['x'] - tx)**2 + (cand['y'] - ty)**2)
            if dist < min_dist:
                min_dist = dist
                best_idx = i
        selected_indices.append(best_idx)
        
    res = env.evaluate(selected_indices, n)
    return res

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    stations_data = os.path.join(current_dir, 'data', 'coverage_results.csv')
    rtree_index = '/home/zoe/Documents/code/Q-view/mergeindex/indexes/ABS_selection.3idx'

    try:
        env = CoverageEnv(stations_data, rtree_path=rtree_index, num_users=600)
        
        print("\n" + "="*70)
        print("      UNIFORM GEOMETRIC PLACEMENT COMPARISON (5-10 Stations)")
        print("="*70)
        print(f"{'N':<4} | {'热点区域覆盖率':<10} | {' 总覆盖率':<10} | {'热点区域平均速率':<10} | {'热点区域最低速率':<10}| {'单机平均负载系数':<10}| {'运行时间：':<10}")
        print("-" * 70)

        results = []
        for n in range(5, 11):
            t1 = time.time()
            res = run_uniform_placement(env, n)
            t2 = time.time()
            print(f"{n:<4} | {res['hs_ratio']*100:>8.2f}% | {res['total_ratio']*100:>8.2f}% | {res['avg_rate']/1e6:>8.3f}Mbps |{res['min_rate'] /1e6:.3f}Mbps |{res['avg_load']*100:.2f}%|{t2-t1}s")
            results.append(res)
            
        print("="*70)
        
    except Exception as e:
        print(f"运行时出错: {e}")
