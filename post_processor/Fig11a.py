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

    return filtered_sorted_data

def smooth_data(data, window_size=40):
    """
    Apply a rolling window smooth for 'f_frac' and 'f_grad' with specified window size.
    Gaussian smoothing is also applied for further smoothing effect.
    """
    data['smooth'] = data['crack_length_derivative'].rolling(window=window_size, min_periods=1, center=True).mean()
    data['smooth'] = gaussian_filter1d(data['smooth'], sigma=5)
    
    return data

def plot_fig11a(data1, data2):
    flex1=load_and_deduplicate_data(data1)
    flex2=load_and_deduplicate_data(data2)
    f0_1 =load_and_deduplicate_data('f0_1.txt')

    # 启用LaTeX渲染和设置全局字体大小
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 24
    plt.rcParams['font.family'] = 'Times New Roman'
    
    fig, ax = plt.subplots(figsize=(10, 8))
    filtered_data = flex1[flex1['crack_length'] <= 90]
    ax.plot(filtered_data['t'], filtered_data['crack_length'], label=r'$P\downarrow$', marker='o', markersize=12, markevery=20, lw=3, color='#008a9d')
    filtered_data = flex2[flex2['crack_length'] <= 90]
    ax.plot(filtered_data['t'], filtered_data['crack_length'], label=r'$P\uparrow$', marker='s', markersize=10, markevery=15, lw=3, color='#eeb0b0')
    ax.plot(f0_1['t'], f0_1['crack_length'], label=r'$P\rightarrow$', ls="--", lw=2, color='black')

    ax.set_ylabel(r'Crack extension $\Delta x$')
    ax.set_xlabel(r'Time step $t$')

    # 设置图例
    ax.legend(loc='best')

    plt.tight_layout()
    plt.savefig('Fig11a.svg')
    plt.show()

    return

plot_fig11a('f0_3.txt', 'f0_4.txt')




