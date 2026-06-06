import time
from scipy.spatial import KDTree
import struct
import numpy as np
import pandas as pd
import os

# 预计算比特计数表（加速bitmap的1的个数统计）
_pop_range = np.arange(256, dtype=np.uint8)
POPCOUNT_TABLE = (_pop_range[:, None] & (1 << np.arange(8)) > 0).sum(1).astype(np.uint8)

class EnvironmentManager:
    def __init__(self, csv_path: str, bin_path: str):
        self.csv_path = csv_path
        self.bin_path = bin_path
        self.station_coords = None
        self.bitmaps = {}        # 基站ID -> 覆盖位图（uint8数组，125000字节）
        self.individual_counts = {} # 基站ID -> 单个基站覆盖的像素数
        self.kd_tree = None
        self.all_station_ids = []
        self.BITMAP_SIZE = 125000  # 1000*1000 / 8 = 125000

    def load_data(self):
        """加载基站坐标和覆盖位图数据"""
        print("[-] Loading spatial data from CSV...")
        try:
            df = pd.read_csv(self.csv_path)
            self.station_coords = df.set_index('station_id')[['x', 'y', 'z']]
            self.all_station_ids = df['station_id'].tolist()
            
            print("[-] Building Spatial Index (KDTree)...")
            self.kd_tree = KDTree(self.station_coords.values)
        except Exception as e:
            print(f"[Error] Failed to load CSV: {e}")
            return

        print("[-] Loading bitmaps from BIN (High Performance Mode)...")
        self._load_binary_bitmaps()

    def _load_binary_bitmaps(self):
        """从二进制文件加载覆盖位图"""
        count = 0
        try:
            with open(self.bin_path, 'rb') as f:
                while True:
                    # 读取基站ID（4字节小端）
                    id_bytes = f.read(4)
                    if not id_bytes: break
                    station_id = struct.unpack('<i', id_bytes)[0]
                    
                    # 读取125KB的位图数据
                    bitmap_bytes = f.read(self.BITMAP_SIZE)
                    if len(bitmap_bytes) != self.BITMAP_SIZE: break
                    
                    if station_id in self.station_coords.index:
                        arr = np.frombuffer(bitmap_bytes, dtype=np.uint8).copy()
                        self.bitmaps[station_id] = arr         
                        # 统计单个基站的覆盖像素数
                        self.individual_counts[station_id] = POPCOUNT_TABLE[arr].sum()
                        
                        count += 1
                        if count % 500 == 0:
                            print(f"    Loaded {count} stations...")
                            
        except Exception as e:
            print(f"[Error] {e}")
            
        print(f"[-] Loaded {len(self.bitmaps)} station bitmaps.")

    def get_neighbors(self, station_id, radius):
        """通过KDTree获取指定半径内的邻居基站ID"""
        if station_id not in self.station_coords.index: return []
        coords = self.station_coords.loc[station_id].values
        indices = self.kd_tree.query_ball_point(coords, r=radius)
        neighbor_ids = [self.all_station_ids[i] for i in indices]
        if station_id in neighbor_ids: neighbor_ids.remove(station_id)
        return neighbor_ids

class GreedySelector:
    def __init__(self, env: EnvironmentManager, n_stations: int):
        self.env = env
        self.n_stations = n_stations  # 需要选择的基站数量
        self.selected_ids = []        # 选中的基站ID列表
        self.total_coverage = 0       # 总覆盖像素数
        self.overlap = 0.0            # 覆盖重叠率
        self.run_time = 0.0           # 算法运行时间

    def calculate_additional_coverage(self, candidate_id, current_union):
        """计算候选基站能新增的覆盖像素数"""
        candidate_bitmap = self.env.bitmaps[candidate_id]
        # 计算新增覆盖：候选位图 - 已覆盖位图（按位与取反后再按位与）
        additional_bitmap = np.bitwise_and(candidate_bitmap, np.bitwise_not(current_union))
        additional_count = POPCOUNT_TABLE[additional_bitmap].sum()
        return additional_count

    def run(self):
        """执行贪心选择算法"""
        print(f"\n[*] Starting Greedy Algorithm for {self.n_stations} stations...")
        start_time = time.time()
        
        # 初始化：已覆盖的位图（全0）
        current_union = np.zeros(self.env.BITMAP_SIZE, dtype=np.uint8)
        remaining_ids = list(self.env.bitmaps.keys())  # 未选中的基站ID
        
        for step in range(self.n_stations):
            if not remaining_ids:
                print(f"[Warning] No more stations to select (only {step} selected)")
                break
            
            best_id = None
            max_additional = -1
            
            # 遍历所有未选中的基站，找到新增覆盖最大的
            for candidate_id in remaining_ids:
                additional = self.calculate_additional_coverage(candidate_id, current_union)
                if additional > max_additional:
                    max_additional = additional
                    best_id = candidate_id
            
            # 选中最优基站
            self.selected_ids.append(best_id)
            remaining_ids.remove(best_id)
            
            # 更新已覆盖位图和总覆盖数
            current_union = np.bitwise_or(current_union, self.env.bitmaps[best_id])
            self.total_coverage = POPCOUNT_TABLE[current_union].sum()
            
            # 计算重叠率：(各基站单独覆盖总和 - 总覆盖) / 总覆盖
            sum_individual = sum(self.env.individual_counts[sid] for sid in self.selected_ids)
            self.overlap = (sum_individual - self.total_coverage) / self.total_coverage if self.total_coverage > 0 else 0.0
            
            # 打印每一步的结果
            print(f"Step {step+1}: Selected ID={best_id}, Additional Coverage={max_additional}, Total Coverage={self.total_coverage}, Overlap={self.overlap:.4f}")
        
        # 计算运行时间
        self.run_time = time.time() - start_time
        
        # 输出最终结果
        print("\n[=] Greedy Optimization Finished.")
        print(f"Final Total Coverage: {self.total_coverage} pixels")
        print(f"Final Overlap Rate: {self.overlap:.4f}")
        print(f"Selected Station IDs: {self.selected_ids}")
        print(f"Greedy Algorithm Run Time: {self.run_time:.6f} seconds")

if __name__ == "__main__":
    # 检查数据文件是否存在
    csv_path = 'data/coverage_results.csv'
    bin_path = 'data/coverage.bin'
    if not os.path.exists(csv_path) or not os.path.exists(bin_path):
        print("[!] Please provide 'coverage_results.csv' and 'coverage.bin' in 'data' folder")
        exit()

    # 1. 加载数据
    load_start = time.time()
    env = EnvironmentManager(csv_path, bin_path)
    env.load_data()
    load_time = time.time() - load_start
    print(f"\nData Loading Time: {load_time:.6f} seconds")

    # 2. 配置贪心算法参数（与GeoAIS对齐）
    n_stations = 60  # 需要选择的无人机基站数量（可修改为其他值，如5、10等）
    
    # 3. 运行贪心算法
    selector = GreedySelector(env, n_stations)
    selector.run()

    # 4. 输出汇总信息（便于对比）
    print("\n===== Greedy Algorithm Summary =====")
    print(f"Number of Selected Stations: {n_stations}")
    print(f"Total Coverage: {selector.total_coverage} pixels")
    print(f"Overlap Rate: {selector.overlap:.4f}")
    print(f"Selected IDs: {selector.selected_ids}")
    print(f"Algorithm Run Time: {selector.run_time:.6f} s")
    print(f"Total Time (Load + Run): {load_time + selector.run_time:.6f} s")