# 优化后的代码
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

def fig11b(filename):
    # 调用函数处理数据和绘制图表
    deduplicated_data = load_and_deduplicate_data(filename)
    intersections = find_intersections(deduplicated_data['t'].values, deduplicated_data['g_frac'].values, deduplicated_data['g_v'].values)
    print("Intersections at time steps:", intersections)
    smoothed_data = smooth_data(deduplicated_data, window_size=5)
    time, g_frac = interpolate_data(smoothed_data['t'].values, smoothed_data['g_frac_smooth'].values)
    _, g_v = interpolate_data(smoothed_data['t'].values, smoothed_data['g_v_smooth'].values)

    plot_fig3(time, g_frac, g_v, smoothed_data, intersections)

    return


def load_and_deduplicate_data(filepath):
    """
    Load data from the specified filepath, sort by 'crack_length' and 't', 
    and remove duplicates keeping the first occurrence for each 'crack_length'.
    """
    # Load and sort data
    data = pd.read_csv(filepath, sep="\t", header=None)
    data.columns = ["t", "crack_length", "g_frac", "g_v", "Fy"]
    data = data[data['t'] > 100]
    data_sorted = data.sort_values(by=["crack_length", "t"])

    # Separate data where 'crack_length' = 40.0
    data_40 = data_sorted[data_sorted['crack_length'] == 40.0]

    # Remove duplicates for other 'crack_length' values
    other_data = data_sorted[data_sorted['crack_length'] != 40.0]
    deduplicated_other_data = other_data.drop_duplicates(subset=["crack_length"], keep="first")

    # Combine and sort data
    modified_data = pd.concat([data_40, deduplicated_other_data]).sort_values(by=["crack_length", "t"])
    return modified_data

    return deduplicated_data

def smooth_data(data, window_size=5):
    """
    Apply a rolling window smooth for 'g_frac' and 'g_v' with specified window size.
    Gaussian smoothing is also applied for further smoothing effect.
    """
    data['g_frac_smooth'] = data['g_frac'].rolling(window=window_size, min_periods=1, center=True).mean()
    data['g_frac_smooth'] = gaussian_filter1d(data['g_frac_smooth'], sigma=5)
    
    data['g_v_smooth'] = data['g_v'].rolling(window=window_size, min_periods=1, center=True).mean()
    data['g_v_smooth'] = gaussian_filter1d(data['g_v_smooth'], sigma=5)
    
    return data

def find_intersections(x, y1, y2):
    """
    Find intersections between two curves y1 and y2, given the common x.
    """
    dy = np.sign(y1 - y2)
    dy_shifted = np.roll(dy, 1)
    intersection_indices = np.where(dy - dy_shifted != 0)[0]
    return x[intersection_indices[1]]

def interpolate_data(x, y):
    """
    Interpolate missing values in the data.
    """
    # 使用插值创建一个具有更多点的新数据集
    f_interp = interp1d(x, y, kind='linear')
    x_new = np.linspace(x.min(), x.max(), 1000)  # 创建一个新的x数组，具有更多的点
    y_new = f_interp(x_new)  # 对y进行插值

    return x_new, y_new

def plot_fig3(time, g_frac, g_v, data,intersections):
    """
    Plot 'g_frac', 'g_v', and 'Fy' against 't' from the deduplicated data.
    """
    # 启用LaTeX渲染和设置全局字体大小
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 24
    plt.rcParams['font.family'] = 'Times New Roman'

    # 选择在曲线上放置标记的数量
    num_markers = 20

    fig, ax1 = plt.subplots(figsize=(10,8))
    
    # 绘制 g_frac 和 g_v
    ax1.plot(time, abs(g_frac), marker='o', markersize=12, markevery=99, lw=3, label=r'$J_{c}$', color='#008a9d')
    ax1.plot(time, abs(g_v), marker='s', markersize=12, markevery=99, lw=3, label=r'$J$', color='#eeb0b0')
    ax1.set_xlabel(r'Time step $t$ ')
    ax1.set_ylabel(r'Configurational forces $J$')

    # Mark intersections
    ax1.axvline(x=1432, color='grey', linestyle='--')

    # 创建第二个坐标轴用于 crack_length
    # ax2 = ax1.twinx()
    # ax2.plot(data['t'], abs(data['crack_length'])-40, label='$\Delta x$', lw=3, color='#243258', linestyle='-.')
    # ax2.set_ylabel(r'Crack extension $\Delta x$')

    # 设置图例
    # lines, labels = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines + lines2, labels + labels2, loc='upper center')
    ax1.legend(loc='best')
    
    plt.xticks([0, 1000, 1432, 2000, 3000])
    plt.tight_layout()
    plt.savefig('Fig11b.svg')
    plt.show()


fig11b('f0_3.txt')
