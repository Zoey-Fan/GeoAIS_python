import numpy as np
import pandas as pd
import random
import math
import copy
import time
import os
from my_rtree import RTree



# --- 配置参数 ---
AREA_SIZE = 1000
B_TOTAL = 20e6        # 总带宽 20MHz
P_TX_DBM = 30.0       # 发射功率 (30dBm = 1W)
N0_DBM = -90.0        # 噪声功率
R_MIN_TH = 0.5e6      # 速率阈值 0.5Mbps (满足基本通话/低质视频)
K_MAX = 50            # 单无人机最大服务用户量
FC_GHZ = 2.0          # 载波频率 2GHz
X_OFFSET = 835000     # 全局坐标偏移 (对接 RTree)
Y_OFFSET = 815000


class Rect3d:
    def __init__(self, minX=0, minY=0, minZ=0, maxX=0, maxY=0, maxZ=0):
        self.min = [minX, minY, minZ]
        self.max = [maxX, maxY, maxZ]

    def __repr__(self):

        return f"Rect3d(min={self.min}, max={self.max})"

# --- 用户生成 ---
def generate_ground_users(num_total=1200,tree=None):
        """生成包含热点区域的用户分布"""
        users = []
        # 热点 1: 市中心高密度区 250*250附近
        for _ in range(int(num_total * 0.35)):
            r = random.gauss(0, 50)
            theta = random.uniform(0, 2 * math.pi)
            x = int(np.clip(250 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(250 + r * math.sin(theta), 0, AREA_SIZE))
            z=float(tree.get_height_3d(x+X_OFFSET,y+Y_OFFSET)or 0.0)
   
            users.append([x, y, z, True])
    
        # 热点 2: 工业园区扩展区 750*750附近
        for _ in range(int(num_total * 0.35)):
            r = random.gauss(0, 80)
            theta = random.uniform(0, 2 * math.pi)
            x = int(np.clip(750 + r * math.cos(theta), 0, AREA_SIZE))
            y = int(np.clip(750 + r * math.sin(theta), 0, AREA_SIZE))
            z=float(tree.get_height_3d(x+X_OFFSET,y+Y_OFFSET)or 0.0)
            users.append([x, y, z, True])
    
        # 普通散布用户
        while len(users) < num_total:
            x = int(random.uniform(0, AREA_SIZE))
            y = int(random.uniform(0, AREA_SIZE))
            z=float(tree.get_height_3d(x+X_OFFSET,y+Y_OFFSET)or 0.0)
            users.append([x, y, z, False])
    
        return pd.DataFrame(users, columns=['x', 'y', 'z', 'is_hotspot'])


# --- 扩展环境管理 ---
class UserEnvManager:
    def __init__(self, csv_path, rtree_path):
        self.stations = pd.read_csv(csv_path)
        self.tree = RTree()
        if os.path.exists(rtree_path):
            self.tree.load_from_file(rtree_path)
        else:
            print(f"[Warning] RTree index {rtree_path} not found!")

    def get_path_loss(self, dist_3d, station_z, mode='LoS'):
        if mode == 'LoS':
            a, b = -28.0 - 20 * math.log10(FC_GHZ), 22.0
        else:
            a, b = 17.5 - 20 * math.log10(40 * math.pi * FC_GHZ / 3), 46 - 7 * math.log10(station_z)
        return -a + b * math.log10(max(dist_3d, 1e-3))

# --- GeoAIS 核心算法 ---
class GeoAIS_UserCov:
    def __init__(self, env: UserEnvManager, users: pd.DataFrame, params: dict):
        self.env = env
        self.users = users
        self.N_pop = params.get('N_pop', 30)
        self.N_stations = params.get('N_stations', 5)
        self.max_iter = params.get('max_iter', 40)
        self.omega = params.get('omega', 0.8)  # 热点区域覆盖权重
        self.memory_cell = None

    

    def _evaluate(self, antibody_indices: set) -> dict:
        """评估选定基站集合的覆盖性能"""
        selected_abs = self.env.stations.iloc[list(antibody_indices)].copy()
        user_list = self.users.values#1000个用户，xyz，true/false  高密度用户和散户
        
        # 1. 搜集每个用户的可连接基站 (不考虑容量)
        #初始化
        potential_links = [[] for _ in range(len(user_list))]
        #orig_row_idx是基站id，selected_abs是基站四元组
        for idx_in_subset, (orig_row_idx, abs_row) in enumerate(selected_abs.iterrows()):
            abs_pos = np.array([abs_row['x'], abs_row['y'], abs_row['z']])
            #遍历所有用户
            for u_idx, user in enumerate(user_list):
                u_pos = np.array([user[0], user[1], user[2]])
                dist_3d = np.linalg.norm(abs_pos - u_pos)

                mode = None
                if dist_3d <= abs_row['nlos_radius']:
                    rect = Rect3d(X_OFFSET+abs_pos[0], Y_OFFSET+abs_pos[1], abs_pos[2], X_OFFSET+u_pos[0], Y_OFFSET+u_pos[1], u_pos[2]+1.5)
                    mode = 'NLoS' if self.env.tree.intersect3d(rect.min, rect.max) else 'LoS'
                elif dist_3d <= abs_row['los_radius']:
                    rect = Rect3d(X_OFFSET+abs_pos[0], Y_OFFSET+abs_pos[1], abs_pos[2], X_OFFSET+u_pos[0], Y_OFFSET+u_pos[1], u_pos[2]+1.5)
                    if not self.env.tree.intersect3d(rect.min, rect.max): mode = 'LoS'
                
                if mode:
                    pl = self.env.get_path_loss(dist_3d, abs_pos[2], mode)
                    snr = P_TX_DBM - pl - N0_DBM
                    potential_links[u_idx].append({'abs_idx': orig_row_idx, 'snr': snr})

        # 2. 关联阶段 (贪婪分配：用户选 SNR 最好且有余量的)
        abs_loads = {idx: 0 for idx in selected_abs.index}
        associations = [None] * len(user_list)
         # 3. 统计指标计算
        covered_count = 0
        hotspot_covered = 0
        total_hotspot = sum(self.users['is_hotspot'])
        rates = []
        for u_idx, links in enumerate(potential_links):
            sorted_links = sorted(links, key=lambda x: x['snr'], reverse=True)
            for link in sorted_links:
                if abs_loads[link['abs_idx']] < K_MAX:
                    abs_loads[link['abs_idx']] += 1
                    bw_user = B_TOTAL / (len(antibody_indices) * abs_loads[link['abs_idx']])
                    snr_lin = 10**(link['snr']/10)
                    rate = bw_user * math.log2(1 + snr_lin)#通信速率
                
                    if rate >= R_MIN_TH:
                        covered_count += 1
                        rates.append(rate)
                        if user_list[u_idx][3]:
                            hotspot_covered += 1
                    break
        
                
        #700个热点用户，300个散户
        hs_ratio = hotspot_covered / total_hotspot if total_hotspot > 0 else 0
        total_ratio = covered_count / len(user_list)
        affinity = self.omega * hs_ratio + (1 - self.omega) * total_ratio
        
        return {
            "antibody": antibody_indices,
            "affinity": affinity,
            "hs_ratio": hs_ratio,
            "total_ratio":total_ratio,
            "avg_load": np.mean(list(abs_loads.values())) / K_MAX,
            "avg_rate": np.mean(rates) if rates else 0,
            "min_rate": np.min(rates) if rates else 0
        }
    
                
    def run(self):
        print(f"[*] Starting User-Centric AIS Optimization...")
        valid_ids = self.env.stations.index.tolist() # 所有可选基站ID 3498个
        #30个个体，每个个体一组基站选择方案（无人机数量）
        population = [set(random.sample(valid_ids, self.N_stations)) for _ in range(self.N_pop)]
    
        for gen in range(self.max_iter):
            eval_pop = [self._evaluate(ind) for ind in population]#调用evaluate给每个个体打分
            eval_pop.sort(key=lambda x: x['affinity'], reverse=True)#按分数从高到低排序
        
            if not self.memory_cell or eval_pop[0]['affinity'] > self.memory_cell['affinity']:
                self.memory_cell = copy.deepcopy(eval_pop[0])
        
            print(f"Gen {gen:2d} | Affinity: {eval_pop[0]['affinity']:.4f} | HotspotCov: {eval_pop[0]['hs_ratio']*100:5.2f}%")
        
            new_pop = []
            # 取前一半优秀个体
            elites = eval_pop[:self.N_pop // 2]
            for rank, elite in enumerate(elites):
                parent = list(elite['antibody'])
            # 关键修复1：确保parent长度等于N_stations（防止集合去重导致长度不足）
                if len(parent) < self.N_stations:
                    # 补充随机基站ID至指定长度
                    supplement = random.sample([x for x in valid_ids if x not in parent], self.N_stations - len(parent))
                    parent += supplement
                elif len(parent) > self.N_stations:
                    # 截断至指定长度
                    parent = parent[:self.N_stations]
                
                num_clones = max(2, (self.N_pop // (rank + 1)))#计算这个优秀方案要复制多少份（排名越靠前，克隆越多））
                for _ in range(num_clones):
                    child = parent.copy()
                    # 关键修复2：基于child实际长度生成随机索引
                    mutate_times = random.randint(1, 2)
                    for _m in range(mutate_times):
                        # 确保索引在child长度范围内
                        idx = random.randint(0, len(child)-1)
                        # 确保新ID不等于原ID（可选优化）
                        new_id = random.choice([x for x in valid_ids if x != child[idx]])
                        child[idx] = new_id
                    # 最后再转 set 去重后，重新补全至N_stations长度
                    child_set = set(child)
                    if len(child_set) < self.N_stations:
                        supplement = random.sample([x for x in valid_ids if x not in child_set], self.N_stations - len(child_set))
                        child_set.update(supplement)
                    new_pop.append(child_set)# 把变异后的新方案加入新一代种群
                    
            while len(new_pop) < self.N_pop:#不够的个体用随机新方案补齐
                new_pop.append(set(random.sample(valid_ids, self.N_stations)))
            population = new_pop[:self.N_pop]#形成新一代种群
        
        return self.memory_cell




if __name__ == "__main__":
    stations_data = './data/coverage_results.csv'
    rtree_index = '/home/zoe/Documents/code/Q-view/mergeindex/indexes/ABS_selection.3idx'
    
    if not os.path.exists(stations_data):
        print("错误: 找不到基站候选集文件。")
        exit()

    env_mod = UserEnvManager(stations_data, rtree_index)
    users_data = generate_ground_users(num_total=600,tree=env_mod.tree)
    #保存到本地
    df=pd.DataFrame(users_data)
    df.to_csv('data/user_csv_GeoAIS.csv',index=False,encoding='utf-8')
    print(f"用户坐标已保存到: {'data/user_csv_GeoAIS.csv'}")
    params_opt = {
        'N_pop': 30,
        'N_stations': 5,
        'max_iter': 200,
        'omega': 0.85
    }
    
    ais_opt = GeoAIS_UserCov(env_mod, users_data, params_opt)
    t1=time.time()
    res = ais_opt.run()
    t2=time.time()
    print("\n" + "="*50)
    print("      MULTI-UAV DEPLOYMENT RESULTS (UserCov)")
    print("="*50)
    print(f"ABS 部署配置 (ID集): {res['antibody']}")
    print(f"热点区域覆盖率:       {res['hs_ratio'] * 100:.2f}%")
    print(f"总覆盖率:            {res['total_ratio']*100:.2f}%")
    print(f"热点区域平均速率:     {res['avg_rate'] / 1e6:.3f} Mbps")
    print(f"热点区域最低速率:     {res['min_rate'] / 1e6:.3f} Mbps")
    print(f"单机平均负载系数:     {res['avg_load'] * 100:.2f}%")
    print(f"运行时间：           {t2-t1}s")
    print("="*50)
