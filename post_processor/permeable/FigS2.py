# 优化后的代码
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

def fig8(data1,data2):
    # 调用函数处理数据和绘制图表
    flex1 = load_and_deduplicate_data(data1)
    flex2 = load_and_deduplicate_data(data2)

    # 添加一个新的 x 坐标列，为索引除以 2
    flex1['x'] = flex1.index / 2
    flex2['x'] = flex2.index / 2
    
    plot_fig8(flex1, flex2)

    return

def load_and_deduplicate_data(filepath):
    """
    Load data from the specified filepath, sort by 'crack_length' and 't', 
    and remove duplicates keeping the first occurrence for each 'crack_length'.
    """
    # Load and sort data
    data = pd.read_csv(filepath, sep="\t", header=None)
    data.columns = ["crack_length"]
    data = data[data['crack_length'] <= 90]
    data_sorted = data.sort_values(by=["crack_length"])

    # Separate data where 'crack_length' = 40.0
    data_40 = data_sorted[data_sorted['crack_length'] == 40.0]

    # Remove duplicates for other 'crack_length' values
    other_data = data_sorted[data_sorted['crack_length'] != 40.0]
    deduplicated_other_data = other_data.drop_duplicates(subset=["crack_length"], keep="last")

    # Combine and sort data
    modified_data = pd.concat([data_40, deduplicated_other_data]).sort_values(by=["crack_length"])
    return deduplicated_other_data

    return deduplicated_data

def plot_fig8(crack_flex1, crack_flex2):
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
    
    ax1.plot(crack_flex1['x'].iloc[::2],crack_flex1['crack_length'].iloc[::2], marker='s', markersize=12, markevery=6, lw=3, label=r'$P\uparrow$', color='#008a9d')
    ax1.plot(crack_flex2['x'].iloc[::2],crack_flex2['crack_length'].iloc[::2], marker='s', markersize=12, markevery=5, lw=3, label=r'$P\downarrow$', color='#eeb0b0')
    ax1.set_xlabel(r'Time step $t$ ')
    ax1.set_ylabel(r'Crack extension $\Delta x$')

    ax1.legend(loc='best')
    # ax1.set_ylim(0,200)

    plt.tight_layout()
    plt.savefig('FigS2.pdf')
    plt.show()

fig8('flex3','flex4')
