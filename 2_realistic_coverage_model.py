import time
import numpy as np
import csv
import os
import struct
import multiprocessing
from my_rtree import RTree


EPSILON = 1e-6
X_OFFSET = 835000
Y_OFFSET = 815000
GRID_WIDTH = 1000
GRID_HEIGHT = 1000
TOTAL_PIXELS = GRID_WIDTH * GRID_HEIGHT
BITMAP_SIZE = TOTAL_PIXELS // 8  # 125,000 bytes

def setup_worker(r_tree):
    """Inital index object"""
    global worker_tree
    worker_tree = r_tree

def generate_discrete_circle(center_x, center_y, radius, resolution=1.0):
    """Generate discrete points inside a circle"""
    min_x = int(np.floor(center_x - radius))
    max_x = int(np.ceil(center_x + radius))
    min_y = int(np.floor(center_y - radius))
    max_y = int(np.ceil(center_y + radius))

    for x in range(min_x, max_x + 1, int(resolution)):
        for y in range(min_y, max_y + 1, int(resolution)):
            distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
            if distance <= radius:
                yield (x, y)

class Rect3d:
    def __init__(self, minX, minY, minZ, maxX, maxY, maxZ):
        self.min = [minX, minY, minZ]
        self.max = [maxX, maxY, maxZ]

def viewshed_analysis_yield(tree, station_id, position, los_points_gen, nlos_points):
    """
    Perform visibility analysis and return coverage point list
    """
    final_covered_points = set()
    nlos_points_set = {tuple(p) for p in nlos_points}

    # NLoS zone
    for point in nlos_points_set:
        #px, py = point
        #pz = tree.get_height_3d(px + X_OFFSET, py + Y_OFFSET)
        #if pz < 0: continue

        #search_rect = Rect3d(position[0] + X_OFFSET, position[1] + Y_OFFSET, position[2], 
        #                     px + X_OFFSET, py + Y_OFFSET, pz)
        #is_hit, _ = tree.intersect3d(search_rect.min, search_rect.max)
         
        final_covered_points.add(point)
        
          

    # LoS zone
    for point_tuple in los_points_gen:
        point = tuple(point_tuple)
        if point in nlos_points_set: continue

        px, py = point
        pz = tree.get_height_3d(px + X_OFFSET, py + Y_OFFSET)
        if pz < 0: continue

        search_rect = Rect3d(position[0] + X_OFFSET, position[1] + Y_OFFSET, position[2], 
                             px + X_OFFSET, py + Y_OFFSET, pz)
        is_hit, _ = tree.intersect3d(search_rect.min, search_rect.max)
        if not is_hit:
            final_covered_points.add(point)

    return list(final_covered_points)

def process_station_task(row):
    """Task processing function for a single base station"""
    try:
        station_id = int(row[0])
        x, y, z = float(row[1]), float(row[2]), float(row[3])
        los_radius, nlos_radius = float(row[4]), float(row[5])

        los_circle_points = generate_discrete_circle(x, y, los_radius)
        nlos_circle_points = generate_discrete_circle(x, y, nlos_radius)

        covered_points = viewshed_analysis_yield(
            worker_tree, station_id, (x, y, z), los_circle_points, nlos_circle_points
        )

        return {"station_id": station_id, "points": covered_points}
    except Exception as e:
        print(f"!! something is wrong with ABS {row[0]}  {e}")
        return None

def point_to_index(x, y):
    """Convert coordinates to bitmap index"""
    ix, iy = int(x), int(y)
    if 0 <= ix < GRID_WIDTH and 0 <= iy < GRID_HEIGHT:
        return iy * GRID_WIDTH + ix
    return None

if __name__ == "__main__":
    # Path definition 
    csv_file = "./data/coverage_results.csv"
    output_bin_file = "./data/coverage.bin"
    index_file = "./data/ABS_selection.3idx"

    # load index
    print("--- Loading Q-view index... ---")
    tree = RTree()
    tree.load_from_file(index_file)

    # read all missions
    tasks_to_run = []
    with open(csv_file, 'r', newline='', encoding='utf-8') as f_in:
        reader = csv.reader(f_in)
        next(reader) 
        for row in reader:
            tasks_to_run.append(row)
    
    total_tasks = len(tasks_to_run)
    print(f"--- successfully read missions:{total_tasks} ---")

    # Parallel computing and binary writing
    num_processes = 8 
    print(f"--- use process number: {num_processes} ... ---")
    
    start_time = time.time()
    processed_count = 0

    with multiprocessing.Pool(processes=num_processes, 
                             initializer=setup_worker, 
                             initargs=(tree,), 
                             maxtasksperchild=1) as pool:
        
        with open(output_bin_file, 'wb') as f_out:
           
            for result in pool.imap_unordered(process_station_task, tasks_to_run):
                if result:
                    station_id = result['station_id']
                    points = result['points']
                    f_out.write(struct.pack('<i', station_id))
                    # generate bitmap
                    bitmap = bytearray(BITMAP_SIZE)
                    for px, py in points:
                        idx = point_to_index(px, py)
                        if idx is not None:
                            byte_idx = idx // 8
                            bit_offset = idx % 8
                            bitmap[byte_idx] |= (1 << bit_offset)

            
                    f_out.write(bitmap)
                    
                    processed_count += 1
                    if processed_count % 10 == 0 or processed_count == total_tasks:
                        print(f"process: {processed_count}/{total_tasks} | ABS {station_id} is already writted")

    end_time = time.time()
    print(f"\n--- Calculate Over！---")
    print(f"Total time: {end_time - start_time:.2f} s")
    print(f"Output file: {output_bin_file} ({os.path.getsize(output_bin_file) / (1024*1024):.2f} MB)")