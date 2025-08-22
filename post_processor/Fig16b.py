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
    data = pd.read_csv(data)
    
    # Filtering data for crack_length != 40 and keeping the first occurrence of each crack_length
    filtered_sorted_data = data.drop_duplicates(subset=['length'], keep='last')

    # Sorting based on 't' to maintain the order
    filtered_sorted_data = filtered_sorted_data.sort_values(by='Time')

    return filtered_sorted_data

def smooth_data(data, window_size=5):
    """
    Apply a rolling window smooth for 'f_frac' and 'f_grad' with specified window size.
    Gaussian smoothing is also applied for further smoothing effect.
    """
    data['smooth'] = data['f_frac'].rolling(window=window_size, min_periods=1, center=True).mean()
    data['smooth'] = gaussian_filter1d(data['smooth'], sigma=5)
    
    return data

def plot_fig16b(data1, data2):
    flex1=load_and_deduplicate_data(data1)
    flex2=load_and_deduplicate_data(data2)

    smooth_flex1=smooth_data(flex1)
    smooth_flex2=smooth_data(flex2)

    # 启用LaTeX渲染和设置全局字体大小
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 24
    plt.rcParams['font.family'] = 'Times New Roman'
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(flex1['length']-40, flex1['smooth'], label=r'$f_{ij}=0$', marker='o', markersize=12, markevery=20, lw=3, color='#008a9d')
    ax.plot(flex2['length']-40, flex2['smooth']-180, label=r'$f_{ij}\neq 0$', marker='s', markersize=10, markevery=15, lw=3, color='#eeb0b0')

    ax.set_xlabel(r'Crack extension $\Delta x$')
    ax.set_ylabel(r'Fracture energy $\psi_{\mathrm{frac}}$')

    # 设置图例
    ax.legend(loc='best')

    plt.tight_layout()
    plt.savefig('Fig16b.svg')
    plt.show()

    return

plot_fig16b('f03_Energy.csv', 'flex3_Energy.csv')




