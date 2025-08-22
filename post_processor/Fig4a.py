import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

def load_and_deduplicate_data(data):
    """
    Load data from the specified filepath, sort by 'crack_length' and 't', 
    and remove duplicates keeping the first occurrence for each 'crack_length'.
    """
    # Load and sort configurational force
    data = pd.read_csv(data, sep="\t", header=None)
    data.columns = ["t", "crack_length", "g_frac", "g_v", "Fy"]
    
    # Filtering data for crack_length != 40 and keeping the first occurrence of each crack_length
    filtered_sorted_data = data.drop_duplicates(subset=['crack_length'], keep='first')

    # Sorting based on 't' to maintain the order
    filtered_sorted_data = filtered_sorted_data.sort_values(by='t')

    filtered_sorted_data.to_csv('f0_1.csv')

    return filtered_sorted_data

def derivative(data):
    #  Compute the differences
    data['time_diff'] = data['t'].diff()
    data['crack_length_diff'] = data['crack_length'].diff()
    # Calculate the derivative
    data['crack_length_derivative'] = data['crack_length_diff'] / data['time_diff']

    # Cleaning up by removing the first row (as its derivative cannot be computed) or filling it with zero or another value
    data = data.dropna().reset_index(drop=True)

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

def smooth_data(data, window_size=40):
    """
    Apply a rolling window smooth for 'f_frac' and 'f_grad' with specified window size.
    Gaussian smoothing is also applied for further smoothing effect.
    """
    data['smooth'] = data['crack_length_derivative'].rolling(window=window_size, min_periods=1, center=True).mean()
    data['smooth'] = gaussian_filter1d(data['smooth'], sigma=5)
    
    return data


def plot_fig4(data1, data2):
    f0_1=load_and_deduplicate_data(data1)
    f0_2=load_and_deduplicate_data(data2)

    df0_1=derivative(f0_1)
    df0_2=derivative(f0_2)

    time, interp_f01_crack_length = interpolate_data(df0_1['t'], df0_1['crack_length'])
    _, interp_f02_crack_length = interpolate_data(df0_2['t'], df0_2['crack_length'])
    f0_1_new=pd.DataFrame({'t':time,'crack_length':interp_f01_crack_length})
    f0_2_new=pd.DataFrame({'t':time,'crack_length':interp_f02_crack_length})
    f0_1_d=derivative(f0_1_new)
    f0_2_d=derivative(f0_2_new)
    f0_1_d=smooth_data(f0_1_d)
    f0_2_d=smooth_data(f0_2_d)

    selected=f0_1_d[0:889]
    coeffs = np.polyfit(selected['t'],selected['crack_length_derivative'], 3)
    p = np.poly1d(coeffs)

    # 启用LaTeX渲染和设置全局字体大小
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 24
    plt.rcParams['font.family'] = 'Times New Roman'
    
    fig, ax1 = plt.subplots(figsize=(10, 8))
    ax1.plot([0,379], [0,0], ls='--', marker='X', markersize=10, lw=3,  color='#b1c44e')
    ax1.plot(f0_2_new['t'], f0_2_new['crack_length']-40, ls='--', marker='X', markersize=10, markevery=99, lw=3, label=r'$P\leftarrow$', color='#b1c44e')
    ax1.plot([0,379], [0,0], ls=':', marker='p', markersize=10, lw=3,  color='#631f66')
    ax1.plot(f0_1_new['t'], f0_1_new['crack_length']-40, ls=':', marker='p', markersize=10, markevery=89, lw=3, label=r'$P\rightarrow$', color='#631f66')
    ax1.set_xlabel(r'Time step $t$')
    ax1.set_ylabel(r'Crack extension $\Delta x$')

    # Mark intersections
    ax1.axvline(x=379, color='grey', linestyle='--')

    # 创建第二个坐标轴用于 crack_length
    ax2 = ax1.twinx()
    ax2.plot(f0_2_d['t'], f0_2_d['smooth'], label=r'$P\leftarrow$', marker='o', markersize=12, markevery=99, lw=3, color='#008a9d')
    ax2.plot(f0_1_d['t'], f0_1_d['smooth'], label=r'$P\rightarrow$', marker='s', markersize=10, markevery=89, lw=3, color='#eeb0b0')

    ax2.set_ylabel(r'Crack extension rate ${\mathrm{d}\Delta x}/{\mathrm{d}t}$')

    # 设置图例
    # lines, labels = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(loc='best',title=r'$\Delta x$')
    ax2.legend(loc='best',title=r'${\mathrm{d}\Delta x}/{\mathrm{d}t}$')

    plt.xticks([0, 379, 1000, 1500, 2000])

    plt.tight_layout()
    plt.savefig('Fig4a.svg')
    plt.show()

    return

plot_fig4('f0_1.txt', 'f0_2.txt')




