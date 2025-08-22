# 优化后的代码
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

def fig16a(data1,data2):
    # 调用函数处理数据和绘制图表
    flex1 = load_and_deduplicate_data(data1)
    flex2 = load_and_deduplicate_data(data2)

    smoothed_flex1 = smooth_data(flex1, window_size=5)
    smoothed_flex2 = smooth_data(flex2, window_size=5)

    time, crack_flex1 = interpolate_data(smoothed_flex1['t'].values, smoothed_flex1['crack_length'].values)
    _, g_v_flex1 = interpolate_data(smoothed_flex1['t'].values, smoothed_flex1['g_v_smooth'].values)

    time, crack_flex2 = interpolate_data(smoothed_flex2['t'].values, smoothed_flex2['crack_length'].values)
    _, g_v_flex2 = interpolate_data(smoothed_flex2['t'].values, smoothed_flex2['g_v_smooth'].values)

    plot_fig16a(time, crack_flex1, crack_flex2, g_v_flex1, g_v_flex2)

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
    data = data[data['crack_length'] <= 90]
    data_sorted = data.sort_values(by=["crack_length", "t"])

    # Separate data where 'crack_length' = 40.0
    data_40 = data_sorted[data_sorted['crack_length'] == 40.0]

    # Remove duplicates for other 'crack_length' values
    other_data = data_sorted[data_sorted['crack_length'] != 40.0]
    deduplicated_other_data = other_data.drop_duplicates(subset=["crack_length"], keep="last")

    # Combine and sort data
    modified_data = pd.concat([data_40, deduplicated_other_data]).sort_values(by=["crack_length", "t"])
    return deduplicated_other_data

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

def interpolate_data(x, y):
    """
    Interpolate missing values in the data.
    """
    # 使用插值创建一个具有更多点的新数据集
    f_interp = interp1d(x, y, kind='linear')
    x_new = np.linspace(x.min(), x.max(), 1000)  # 创建一个新的x数组，具有更多的点
    y_new = f_interp(x_new)  # 对y进行插值

    return x_new, y_new

def plot_fig16a(time,crack_flex1, crack_flex2, g_v_flex1, g_v_flex2):
    """
    Plot 'g_frac', 'g_v' against 't' from the deduplicated data.
    """
    # 启用LaTeX渲染和设置全局字体大小
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 24
    plt.rcParams['font.family'] = 'Times New Roman'

    # 选择在曲线上放置标记的数量
    num_markers = 20

    fig, ax1 = plt.subplots(figsize=(10,8))
    
    # 绘制 g_frac 和 g_v
    ax1.plot(crack_flex1-40, abs(g_v_flex1), marker='s', markersize=12, markevery=99, lw=3, label=r'$f_{ij}=0$', color='#008a9d')
    ax1.plot(crack_flex2-40, abs(g_v_flex2), marker='s', markersize=12, markevery=99, lw=3, label=r'$f_{ij}\neq 0$', color='#eeb0b0')
    ax1.set_xlabel(r'Crack extension $\Delta x$')
    ax1.set_ylabel(r'Configurational forces $J$')

    ax1.legend(loc='best')

    plt.tight_layout()
    plt.savefig('Fig16a.svg')
    plt.show()

fig16a('f0_3.txt','flex3.txt')
