import numpy as np
import pandas as pd
import struct
import random
import math
import copy
from scipy.spatial import KDTree
import os
import time


_pop_range = np.arange(256, dtype=np.uint8)
POPCOUNT_TABLE = (_pop_range[:, None] & (1 << np.arange(8)) > 0).sum(1).astype(np.uint8)

# data loading 

class EnvironmentManager:
    def __init__(self, csv_path: str, bin_path: str):
        self.csv_path = csv_path
        self.bin_path = bin_path
        self.station_coords = None
        self.bitmaps = {}        # Dict: station_id -> np.array(uint8) [125000 bytes]
        self.individual_counts = {} # Dict: station_id -> int 
        self.kd_tree = None
        self.all_station_ids = []
        
        self.BITMAP_SIZE = 125000  # 1000*1000 / 8

    def load_data(self):
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
        count = 0
        try:
            with open(self.bin_path, 'rb') as f:
                while True:
                    # read ID
                    id_bytes = f.read(4)
                    if not id_bytes: break
                    station_id = struct.unpack('<i', id_bytes)[0]
                    
                    # ead 125KB bitmap
                    bitmap_bytes = f.read(self.BITMAP_SIZE)
                    if len(bitmap_bytes) != self.BITMAP_SIZE: break
                    
                    if station_id in self.station_coords.index:
                        arr = np.frombuffer(bitmap_bytes, dtype=np.uint8).copy()
                        
                        self.bitmaps[station_id] = arr         
                        self.individual_counts[station_id] = POPCOUNT_TABLE[arr].sum()
                        
                        count += 1
                        if count % 500 == 0:
                            print(f"    Loaded {count} stations...")
                            
        except Exception as e:
            print(f"[Error] {e}")
            
        print(f"[-] Loaded {len(self.bitmaps)} station bitmaps.")

    def get_neighbors(self, station_id, radius):
        """KDTree查询保持不变"""
        if station_id not in self.station_coords.index: return []
        coords = self.station_coords.loc[station_id].values
        indices = self.kd_tree.query_ball_point(coords, r=radius)
        neighbor_ids = [self.all_station_ids[i] for i in indices]
        if station_id in neighbor_ids: neighbor_ids.remove(station_id)
        return neighbor_ids

# 2. GeoAIS
class GeoAIS:
    def __init__(self, env: EnvironmentManager, params: dict):
        self.env = env
        self.N_pop = params.get('N_pop', 50)
        self.N_stations = params.get('N_stations', 20)
        self.max_iter = params.get('max_iter', 100)
        self.beta_clone = params.get('beta_clone', 10)
        self.n_elite = params.get('n_elite', 10)
        self.r_edit = params.get('r_edit', 0.1)
        self.omega = params.get('omega', 0.7)
        self.R_base = params.get('R_base', 1000.0)
        self.gamma = params.get('gamma', 2.0)
        
        self.population = []
        self.memory_cell = None
        self.history = []

    def _init_population(self):
        self.population = []
        valid_ids = list(self.env.bitmaps.keys()) 
        for _ in range(self.N_pop):
            ind = set(random.sample(valid_ids, self.N_stations))
            self.population.append(ind)

    def _evaluate(self, antibody: set) -> dict:
        """
        Innovation 1: Spatial Overlap Penalty
        """
        station_ids = list(antibody)
        arrays = [self.env.bitmaps[sid] for sid in station_ids]
        
        if not arrays:
            return {"antibody": antibody, "f_cov": 0, "phi_overlap": 0, "affinity": 0}
        
        union_bitmap = np.bitwise_or.reduce(arrays)
        
        f_cov = POPCOUNT_TABLE[union_bitmap].sum()
        
        sum_cov = sum(self.env.individual_counts[sid] for sid in station_ids)
        
        if f_cov > 0:
            phi_overlap = (sum_cov - f_cov) / f_cov
        else:
            phi_overlap = 0.0
            
        return {
            "antibody": antibody,
            "f_cov": f_cov,
            "phi_overlap": phi_overlap,
            "affinity": 0.0 
        }

    def _calculate_affinity_batch(self, evaluated_pop: list):
        """
        E(A) = w * Norm(Cov) - (1-w) * Norm(Overlap)
        """
        if not evaluated_pop: return
        
        f_vals = [p['f_cov'] for p in evaluated_pop]
        phi_vals = [p['phi_overlap'] for p in evaluated_pop]
        
        f_max, f_min = max(f_vals), min(f_vals)
        phi_max, phi_min = max(phi_vals), min(phi_vals)
        
        f_range = f_max - f_min if f_max != f_min else 1.0
        phi_range = phi_max - phi_min if phi_max != phi_min else 1.0
        
        for p in evaluated_pop:
            norm_f = (p['f_cov'] - f_min) / f_range
            norm_phi = (p['phi_overlap'] - phi_min) / phi_range
            p['affinity'] = self.omega * norm_f - (1 - self.omega) * norm_phi

    def _geo_hypermutation(self, parent: dict) -> set:
        """
        Innovation 2: Affinity-Modulated Geo-Spatial Hypermutation
        """

        antibody = list(parent['antibody'])
        norm_affinity = max(0, min(1, parent['affinity']))
        r_mut = self.R_base * math.exp(-self.gamma * norm_affinity)# R_mut = R_base * exp(-gamma * E)
        
        new_antibody = set(antibody)
        if not new_antibody: return new_antibody
        station_to_move = random.choice(antibody)  
        neighbors = self.env.get_neighbors(station_to_move, r_mut)       
    
        valid_neighbors = [nid for nid in neighbors if nid in self.env.bitmaps]
        
        if valid_neighbors:
            new_station = random.choice(valid_neighbors)
            if new_station not in new_antibody:
                new_antibody.remove(station_to_move)
                new_antibody.add(new_station)
        
        return new_antibody

    def run(self):
        print(f"[*] Starting GeoAIS (Bitwise Mode) for {self.N_stations} stations...")
        self._init_population()
        self.memory_cell = None
        
        self.cov_history = [] 

        for gen in range(self.max_iter):
            eval_pop = [self._evaluate(ind) for ind in self.population]
            self._calculate_affinity_batch(eval_pop)
            eval_pop.sort(key=lambda x: x['affinity'], reverse=True)
            
            current_best = eval_pop[0]
            
            update_memory = False
            if self.memory_cell is None:
                update_memory = True
            else:
                if current_best['f_cov'] > self.memory_cell['f_cov']:
                    update_memory = True
                elif (current_best['f_cov'] == self.memory_cell['f_cov']) and \
                     (current_best['phi_overlap'] < self.memory_cell['phi_overlap']):
                    update_memory = True
            
            if update_memory:
                self.memory_cell = copy.deepcopy(current_best)
            print(f"Gen {gen}: [New Best] Cov={current_best['f_cov']}, Overlap={current_best['phi_overlap']:.4f},affinity={current_best['affinity']:.6f}")
            

            
            # clone selection
            elites = eval_pop[:self.n_elite]
            clone_pool = []
            for rank, elite in enumerate(elites):
                
                n_clones = int(self.beta_clone * self.N_pop / (rank + 1) / 2)
                n_clones = max(1, n_clones)
                for _ in range(n_clones):
                    mutated = self._geo_hypermutation(elite)
                    clone_pool.append(mutated)
                  
            n_edit = int(self.N_pop * self.r_edit)
            new_blood = []
            valid_ids = list(self.env.bitmaps.keys())
            for _ in range(n_edit):
                new_blood.append(set(random.sample(valid_ids, self.N_stations)))
            
            combined_pop_sets = [p['antibody'] for p in eval_pop] + clone_pool + new_blood
            
            combined_eval = [self._evaluate(ind) for ind in combined_pop_sets]
            self._calculate_affinity_batch(combined_eval) 
            combined_eval.sort(key=lambda x: x['affinity'], reverse=True)
            
            self.population = [p['antibody'] for p in combined_eval[:self.N_pop]]
        
        print("\n[=] Optimization Finished.")
        print(f"Final Solution: Cov={self.memory_cell['f_cov']}, Overlap={self.memory_cell['phi_overlap']:.4f}")
        return self.memory_cell

# 3. main process
if __name__ == "__main__":
    # check
    if not os.path.exists('data/coverage_results.csv') or not os.path.exists('data/coverage.bin'):
        print("[!] Please provide 'coverage_results.csv' and 'coverage.bin'")
        exit()

    # load data
    startTime=time.time()
    env = EnvironmentManager('data/coverage_results.csv', 'data/coverage.bin')
    env.load_data()
    endTime=time.time()
    t=endTime-startTime
    print(f"loading time:{t:.6f} s")
    
    # running
    
    params = {
        'N_pop': 50,      # population
        'N_stations': 3, # number of station needed 
        'max_iter': 100,  # generation
        'beta_clone': 10, # clone
        'n_elite': 10,    # elite number
        'R_base': 1000.0, # the inital radius
        'gamma': 3.0,     
        'omega': 0.7      # coverage rartio
    }
    start_time=time.time()
    algo = GeoAIS(env, params)
    best = algo.run()
    end_time=time.time()
    elapsed_time = end_time - start_time
    print("Final Solution IDs:", best['antibody'])
    print(f"running time: {elapsed_time:.6f} s")
