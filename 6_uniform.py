import numpy as np
import pandas as pd
import math
import os
import struct
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from collections import defaultdict

# ===================== 基础配置 =====================
AREA_SIZE = 1000          # 区域大小 1000x1000m
CSV_FILE_PATH = "./data/coverage_results.csv"  # 基站候选集CSV
BITMAP_FILE_PATH = "./data/coverage.bin"       # 覆盖位图二进制文件
OUTPUT_DIR = "./data/figure"                   # 可视化结果保存目录

# Bitmap解析配置
GRID_WIDTH = 1000
GRID_HEIGHT = 1000
TOTAL_PIXELS = GRID_WIDTH * GRID_HEIGHT
BITMAP_SIZE = TOTAL_PIXELS // 8  # 125,000 bytes

# 要分析的基站数量列表
BS_COUNTS = [5, 10, 20, 30, 40, 50, 60]

# ===================== 工具函数 =====================
def get_uniform_targets(n):
    """生成n个均匀分布的目标坐标（1000x1000区域内）"""
    if n == 5:
        return [(500, 500), (250, 250), (250, 750), (750, 250), (750, 750)]
    
    # 通用网格均匀分布逻辑
    rows = int(math.sqrt(n))
    cols = math.ceil(n / rows)
    
    targets = []
    for r in range(rows):
        for c in range(cols):
            if len(targets) < n:
                x = (c + 0.5) * (AREA_SIZE / cols)
                y = (r + 0.5) * (AREA_SIZE / rows)
                targets.append((x, y))
    return targets

def match_closest_candidate(targets, candidates):
    """根据目标坐标匹配候选集中最近的基站ID"""
    selected_indices = []
    selected_ids = []
    
    for tx, ty in targets:
        best_idx = -1
        min_dist = float('inf')
        best_id = -1
        
        for i, cand in enumerate(candidates):
            if i in selected_indices:
                continue
            dist = math.sqrt((cand['x'] - tx)**2 + (cand['y'] - ty)**2)
            if dist < min_dist:
                min_dist = dist
                best_idx = i
                best_id = cand['station_id']
        
        selected_indices.append(best_idx)
        selected_ids.append(best_id)
    
    return selected_ids, selected_indices

def point_to_index(x, y):
    """坐标转位图索引"""
    ix, iy = int(x), int(y)
    if 0 <= ix < GRID_WIDTH and 0 <= iy < GRID_HEIGHT:
        return iy * GRID_WIDTH + ix
    return None

def index_to_point(idx):
    """位图索引转回坐标"""
    if 0 <= idx < TOTAL_PIXELS:
        iy = idx // GRID_WIDTH
        ix = idx % GRID_WIDTH
        return ix, iy
    return None

def read_station_bitmap(station_id, bin_file_path):
    """从二进制文件读取指定基站的覆盖点集"""
    covered_points = []
    station_id_bytes = struct.pack('<i', station_id)
    
    with open(bin_file_path, 'rb') as f:
        while True:
            id_buf = f.read(4)
            if not id_buf:
                break
            
            bitmap_buf = f.read(BITMAP_SIZE)
            if len(bitmap_buf) != BITMAP_SIZE:
                continue
            
            current_id = struct.unpack('<i', id_buf)[0]
            if current_id == station_id:
                # 解析bitmap
                for byte_idx in range(BITMAP_SIZE):
                    byte = bitmap_buf[byte_idx]
                    if byte == 0:
                        continue
                    for bit_offset in range(8):
                        if byte & (1 << bit_offset):
                            pixel_idx = byte_idx * 8 + bit_offset
                            x, y = index_to_point(pixel_idx)
                            covered_points.append((x, y))
                break
    
    return covered_points

def calculate_coverage_metrics(station_coverage):
    """计算覆盖率和重叠率"""
    point_count = defaultdict(int)
    all_points = []
    
    # 统计每个点被覆盖的次数
    for sid, points in station_coverage.items():
        all_points.extend(points)
        for p in points:
            point_count[p] += 1
    
    # 计算指标
    total_covered = len(point_count)                  # 去重覆盖点数
    total_pixels = GRID_WIDTH * GRID_HEIGHT          # 总像素数
    overlap_points = sum(1 for p, cnt in point_count.items() if cnt >= 2)  # 重叠点数
    
    coverage_rate = total_covered / total_pixels     # 覆盖率
    overlap_rate = overlap_points / total_covered if total_covered > 0 else 0.0  # 重叠率
    
    return {
        'total_covered': total_covered,
        'overlap_points': overlap_points,
        'coverage_rate': coverage_rate,
        'overlap_rate': overlap_rate,
        'point_count': point_count
    }

def visualize_coverage(station_ids, station_info, station_coverage, output_path):
    """可视化基站覆盖范围"""
    # 创建输出目录
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 创建画布
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # 定义颜色映射
    light_colors = [
        mcolors.to_rgba('lightblue', alpha=0.8),
        mcolors.to_rgba('lightgreen', alpha=0.8),
        mcolors.to_rgba('lightcoral', alpha=0.8),
        mcolors.to_rgba('lightyellow', alpha=0.8),
        mcolors.to_rgba('lavender', alpha=0.8),
        mcolors.to_rgba('orange', alpha=0.8)
    ]
    color_map = {sid: light_colors[i % len(light_colors)] for i, sid in enumerate(station_ids)}
    
    # 绘制每个基站的覆盖范围和位置
    for sid in station_ids:
        if sid not in station_info or sid not in station_coverage:
            continue
        
        x, y, z, los_r, nlos_r = station_info[sid]
        coverage_points = station_coverage[sid]
        color = color_map[sid]
        
        # 绘制覆盖点
        if coverage_points:
            cov_x = [p[0] for p in coverage_points]
            cov_y = [p[1] for p in coverage_points]
            ax.scatter(cov_x, cov_y, c=color, alpha=0.6, s=1)
        
        # 绘制基站位置
        ax.scatter(x, y, marker='*', s=200, edgecolors='black', linewidth=1.5)
        
        # 绘制LoS/NLoS半径圆
        los_circle = Circle((x, y), los_r, color=color, alpha=0.2, linestyle='--', linewidth=1)
        nlos_circle = Circle((x, y), nlos_r, color=color, alpha=0.3, linestyle=':', linewidth=1)
        ax.add_patch(los_circle)
        ax.add_patch(nlos_circle)
    
    # 构造图例
    colors = [color_map[sid] for sid in station_ids if sid in color_map]
    legend_points = [Line2D([0], [0], marker='.', color='w', markerfacecolor=colors[i], markersize=9) 
                     for i in range(min(3, len(colors)))]  # 最多显示3个颜色示例
    h_abs = Line2D([0], [0], marker='*', color='w', markerfacecolor='white', 
                   markeredgecolor='black', markersize=12, label='ABS')
    
    handles = [tuple(legend_points), h_abs]
    labels = ['不同ABS覆盖点', 'ABS基站']
    
    # 计算覆盖指标
    metrics = calculate_coverage_metrics(station_coverage)
    
    # 设置图表属性
    ax.set_title(
        f'均匀分布基站覆盖分析（数量：{len(station_ids)}）\n'
        f'覆盖率: {metrics["coverage_rate"]*100:.1f}% | 重叠率: {metrics["overlap_rate"]*100:.1f}%',
        fontsize=14, fontweight='bold'
    )
    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 1000)
    ax.set_xlabel('X (m)', fontsize=12)
    ax.set_ylabel('Y (m)', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    # 添加图例
    ax.legend(
        handles=handles,
        labels=labels,
        loc='lower left',
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.5)},
        handletextpad=0.8,
        fontsize=10
    )
    
    # 保存图片
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"可视化结果已保存至: {output_path}")

def read_station_coords(csv_file_path, target_ids):
    """读取基站坐标和半径信息"""
    df = pd.read_csv(csv_file_path)
    target_df = df[df['station_id'].isin(target_ids)]
    station_info = {}
    
    for _, row in target_df.iterrows():
        sid = int(row['station_id'])
        station_info[sid] = (
            float(row['x']),
            float(row['y']),
            float(row['z']),
            float(row['los_radius']),
            float(row['nlos_radius'])
        )
    
    return station_info

# ===================== 主函数 =====================
def main():
    # 1. 加载基站候选集
    if not os.path.exists(CSV_FILE_PATH):
        raise FileNotFoundError(f"候选集文件不存在: {CSV_FILE_PATH}")
    
    df_candidates = pd.read_csv(CSV_FILE_PATH)
    candidates = df_candidates.to_dict('records')
    print(f"成功加载 {len(candidates)} 个基站候选点")
    
    # 2. 存储所有统计结果
    all_results = []
    
    # 3. 遍历每个基站数量
    for n in BS_COUNTS:
        print(f"\n===== 分析 {n} 个均匀分布基站 =====")
        
        # 生成均匀目标坐标
        targets = get_uniform_targets(n)
        
        # 匹配最近的基站ID
        selected_ids, _ = match_closest_candidate(targets, candidates)
        print(f"匹配到的基站ID: {selected_ids}")
        
        # 读取基站坐标信息
        station_info = read_station_coords(CSV_FILE_PATH, selected_ids)
        
        # 读取每个基站的覆盖点集
        station_coverage = {}
        for sid in selected_ids:
            points = read_station_bitmap(sid, BITMAP_FILE_PATH)
            station_coverage[sid] = points
            print(f"基站 {sid} 覆盖点数: {len(points)}")
        
        # 计算覆盖指标
        metrics = calculate_coverage_metrics(station_coverage)
        print(f"覆盖率: {metrics['coverage_rate']*100:.2f}%")
        print(f"重叠率: {metrics['overlap_rate']*100:.2f}%")
        
        # 保存统计结果
        all_results.append({
            'bs_count': n,
            'station_ids': selected_ids,
            'total_covered': metrics['total_covered'],
            'overlap_points': metrics['overlap_points'],
            'coverage_rate': metrics['coverage_rate'],
            'overlap_rate': metrics['overlap_rate']
        })
        
        # 可视化
        output_path = os.path.join(OUTPUT_DIR, f"uniform_bs_{n}_coverage.png")
        visualize_coverage(selected_ids, station_info, station_coverage, output_path)
    
    # 4. 保存统计结果到CSV
    results_df = pd.DataFrame(all_results)
    results_path = os.path.join(OUTPUT_DIR, "uniform_bs_coverage_stats.csv")
    results_df.to_csv(results_path, index=False, encoding='utf-8')
    print(f"\n所有统计结果已保存至: {results_path}")
    
    # 5. 打印汇总表
    print("\n===== 汇总结果 =====")
    print(f"{'基站数量':<10} {'覆盖率(%)':<12} {'重叠率(%)':<12}")
    print("-" * 35)
    for res in all_results:
        print(f"{res['bs_count']:<10} {res['coverage_rate']*100:<12.2f} {res['overlap_rate']*100:<12.2f}")

if __name__ == "__main__":
    main()