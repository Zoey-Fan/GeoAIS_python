from my_rtree import RTree
import time
from scipy.stats import norm
import numpy as np
import pandas as pd
import math

EPSILON = 1e-6
MAXH = 1e6  # Replace with an appropriate maximum height
MINH = -1e6  # Replace with an appropriate minimum height
tree = None
X=835000
Y=815000

class Rect3d:
    def __init__(self, minX=0, minY=0, minZ=0, maxX=0, maxY=0, maxZ=0):
        self.min = [minX, minY, minZ]
        self.max = [maxX, maxY, maxZ]

    def __repr__(self):

        return f"Rect3d(min={self.min}, max={self.max})"


def generate_station_coordinates(x_range, y_range, z_range):
    stations = []
    x_start, x_end, x_step = x_range
    y_start, y_end, y_step = y_range
    z_start, z_end, z_step = z_range

    for x in range(x_start, x_end + 1, x_step):
        for y in range(y_start, y_end + 1, y_step):
            for z in range(z_start, z_end + 1, z_step):
                stations.append((x, y, z))
    
    return stations

def filter_stations(all_stations):
    
    stations_by_xy = {}
    for x,y,z in all_stations:
        if (x,y) not in stations_by_xy:
            stations_by_xy[(x,y)]=[]
        stations_by_xy[(x,y)].append(z)
    valid_stations=[]
    for (x,y),z_values in stations_by_xy.items():
        TerrainHeight=tree.get_height_3d(X+x,Y+y)
        #print(f"TerrainHeight={TerrainHeight}")
        for z in z_values:
            if z>TerrainHeight and TerrainHeight>0:
                valid_stations.append((x,y,z))

    return valid_stations

def calculate_coverage_radius(h_u,
    link_type,
    path_loss_params,
    snr_threshold_db,
    tx_power_dbm,
    noise_power_dbm,
    coverage_probability):
    if link_type not in path_loss_params:
        raise ValueError("link_type must be LoS or NLoS!")
    params=path_loss_params[link_type]
    a=params['a']
    b=params['b']
    sigma=params['sigma']

    effective_tx_power_db = tx_power_dbm - noise_power_dbm
    gamma_0_bar = snr_threshold_db - effective_tx_power_db
    
    q_inv_epsilon = norm.ppf(1 - coverage_probability)
    numerator = 2*(sigma * q_inv_epsilon - gamma_0_bar + a) * np.log(10)
   
    exponent_term = numerator / b
    
    slant_radius_sq = np.exp(exponent_term)-h_u**2
    slant_radius=math.sqrt(slant_radius_sq)
    return slant_radius

if __name__ == "__main__":
    
    #---S1: Generate the ABS candidate points---
    X_AXIS_RANGE = (0, 1000, 50)
    Y_AXIS_RANGE = (0, 1000, 50)
    Z_AXIS_RANGE = (40, 130, 10)
   
    base_stations_inital = generate_station_coordinates(X_AXIS_RANGE, Y_AXIS_RANGE, Z_AXIS_RANGE)
  
    print(f"successfully generate the ABS candidate points！")
    print(f"The total number is : {len(base_stations_inital)}")
    print("-----------------------------------")  
    
    tree=RTree()
    index="./data/ABS_selection.3idx"
    t1=time.time()
    tree.load_from_file(index)
    t2=time.time()
    print("Index load time:", t2 - t1)
    base_stations=filter_stations(base_stations_inital)
    print(f"successfully filter the ABS candidate points！")
    print(f"The total number is : {len(base_stations)}")
    print("-----------------------------------")  

    #---S2 Calculate the Radius of los and nlos for each base station---
    fc_ghz=2.0 #载波频率 2GHz
    snr_threshold = 25.0         # SNR threshold (dB)
    tx_power = 30.0              # transmission power (dBm)
    noise_power = -90.0          # noise power (dBm)
    cov_prob = 0.9               # coverage ratio threshold
    # path loss parameters (3GPP UMa-AV model)
    path_loss_parameters = {
        'LoS': {
            'a': -28.0 - 20 * np.log10(fc_ghz),
            'b': 22.0,
            'sigma': 4.64
        },
        'NLoS': {
    
            'a': 17.5 - 20 * np.log10(40 * np.pi * fc_ghz / 3),
            'b': 0, # update dynamicly
            'sigma': 6.0
        }
    }
    
    coverage_radius=[]

    for i in range(len(base_stations)):
        x,y,z=base_stations[i]

        path_loss_parameters['NLoS']['b']=46 - 7 * np.log10(z)
        # --- LoS radius ---
        los_slant_radius = calculate_coverage_radius(z,
            link_type='LoS',
            path_loss_params=path_loss_parameters,
            snr_threshold_db=snr_threshold,
            tx_power_dbm=tx_power,
            noise_power_dbm=noise_power,
            coverage_probability=cov_prob
        )

        # --- NLoS radius ---
        nlos_slant_radius = calculate_coverage_radius(z,
            link_type='NLoS',
            path_loss_params=path_loss_parameters,
            snr_threshold_db=snr_threshold,
            tx_power_dbm=tx_power,
            noise_power_dbm=noise_power,
            coverage_probability=cov_prob
        )
    
        coverage_radius.append({
        'station_id': i,
        'x': x,
        'y': y, 
        'z': z,
        'los_radius': los_slant_radius,
        'nlos_radius': nlos_slant_radius
    })
    # 创建DataFrame
    df_results = pd.DataFrame(coverage_radius)
    # 高效存储到文件
    df_results.to_csv('./data/coverage_results.csv', index=False)  
    





