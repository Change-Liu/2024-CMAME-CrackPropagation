import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from statsmodels.nonparametric.smoothers_lowess import lowess

# Set the font to Times New Roman
# Enable LaTeX in Matplotlib
plt.rcParams['text.usetex'] = True
# plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman'] # or another LaTeX compatible font
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24  # You can also set the global font size

# Function to remove outliers using moving median
def remove_outliers_with_rolling(data, window_size=5):
    rolling_median = data.rolling(window=window_size, center=True).median()
    diff = np.abs(data - rolling_median)
    threshold = rolling_median.median() + 1.5 * (rolling_median.quantile(0.75) - rolling_median.quantile(0.25))
    outliers = diff > threshold
    data_no_outliers = data.mask(outliers).interpolate(limit_direction='both')
    return data_no_outliers

# Function to apply LOWESS smoothing
def apply_lowess_smoothing(data, frac=0.1):
    lowess_smoothed = lowess(data['g_diff_no_outliers'], data['crack_length'], frac=frac)
    return lowess_smoothed[:, 1]  # Extract the smoothed y values

# Function to process data, remove outliers and apply smoothing
def process_data(file_path):
    data = pd.read_csv(file_path, delim_whitespace=True, header=None, 
                       names=['time_step', 'crack_length', 'g_frac', 'g_v'])
    data['g_v_no_outliers'] = remove_outliers_with_rolling(data['g_v'])
    data['g_frac_no_outliers'] = remove_outliers_with_rolling(data['g_frac'])
    data['g_diff_no_outliers'] = data['g_v_no_outliers'] - data['g_frac_no_outliers']
    data['lowess_smooth_g_diff'] = apply_lowess_smoothing(data)
    return data

# Function to truncate data starting from the last point where crack length is 40
def truncate_data(data):
    last_index = data[data['crack_length'] == 40].index[-1]
    return data.loc[last_index:].reset_index(drop=True)

# Set font size for all text in the plot
plt.rcParams.update({'font.size': 24})

# File paths to process
files_to_process = ['f0_1.txt', 'f0_2.txt', 'flex1.txt', 'flex2.txt']

# Process and truncate the data for all files
processed_data = {file_path: truncate_data(process_data(file_path)) for file_path in files_to_process}

# Plotting the truncated and smoothed data
plt.figure(figsize=(12, 10)) 
colors = ['c', 'g', 'r', 'b']  # Define the color scheme
labels=[r'$P\leftarrow$', r'$P\rightarrow$', r'with flex $P\leftarrow$', r'with flex $P\rightarrow$']

for i, file_path in enumerate(files_to_process):
    data = processed_data[file_path]
    plt.plot(data['crack_length'], -data['lowess_smooth_g_diff'], label=labels[i],
             color=colors[i], antialiased=True,
             lw=2)

# Configure the plot
plt.xlabel('Crack Length')
plt.ylabel('$g_d^*$')
plt.legend()
plt.gca().spines['top'].set_visible(True)
plt.gca().spines['right'].set_visible(True)
plt.title('')
plt.grid(False)

# Save the figure with truncation and LOWESS smoothing
plt.tight_layout()
plt.savefig('12c-g.svg')

# Show the plot
plt.show()
