import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
import pandas as pd

def load_and_process_csv(csv_file_path):
    df = pd.read_csv(csv_file_path, header=None)
    matrix = np.matrix(df.values)
    matrix = matrix[:, :]
    return matrix

def compute_psd(matrix, fs=2000):
    nperseg = matrix.shape[1]
    psds = []
    print(matrix.shape)
    for row in matrix:
        frequencies, power_density = welch(row, fs=fs, nperseg=nperseg, scaling='density')
        normalized_power_density = power_density / np.sum(power_density)
        psds.append(normalized_power_density)
    
    psds = np.array(psds)

    # Use np.squeeze() to remove unnecessary dimensions
    psds = np.squeeze(psds)
    
    


    
    psd_mean = np.mean(psds, axis=0)
    psd_std = np.std(psds, axis=0)
    #psd_mean = psd_mean/np.max(psd_mean)
    #psd_std = psd_std/np.max(psd_mean)

    return frequencies, psd_mean, psd_std, psds
csv_file_path_1 = '/gpfs/data/doiron-lab/draco/results_nn/exp_8957510/all_raw_activities.csv'
csv_file_path_2 = '/gpfs/data/doiron-lab/draco/results_nn/exp_8957562/all_raw_activities.csv'

# Load, process, and compute PSD for both datasets
matrix_1 = load_and_process_csv(csv_file_path_1)
frequencies_1, psd_mean_1, psd_std_1, psds_1 = compute_psd(matrix_1)

matrix_2 = load_and_process_csv(csv_file_path_2)
frequencies_2, psd_mean_2, psd_std_2, psds_2 = compute_psd(matrix_2)

# Plotting Average Normalized Power Spectral Density
plt.figure(figsize=(10, 6))
plt.plot(frequencies_1[0:10], psd_mean_1[0:10], label='P PSD', color='green')
plt.fill_between(frequencies_1[0:10], psd_mean_1[0:10]-psd_std_1[0:10], psd_mean_1[0:10]+psd_std_1[0:10], color='green', alpha=0.1)
plt.plot(frequencies_2[0:10], psd_mean_2[0:10], label='No P PSD', color='purple')
plt.fill_between(frequencies_2[0:10], psd_mean_2[0:10]-psd_std_2[0:10], psd_mean_2[0:10]+psd_std_2[0:10], color='purple', alpha=0.1)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Normalized Power')
plt.title('Average Normalized Power Spectral Density')
plt.legend()
plt.grid(False)
plt.savefig('combined_average_PSD.png')

# Plotting heat maps
# Dataset 1 Heat Map
plt.figure(figsize=(10, 6))
plt.imshow(psds_1, aspect='auto', cmap='viridis', extent=[frequencies_1[0], frequencies_1[-1], 0, psds_1.shape[0]])
plt.colorbar(label='Normalized Power')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Row in Dataset')
plt.title('Heat Map of Power Spectral Densities for Dataset 1')
plt.savefig('HP_P.png')

# Dataset 2 Heat Map
plt.figure(figsize=(10, 6))
plt.imshow(psds_2, aspect='auto', cmap='viridis', extent=[frequencies_2[0], frequencies_2[-1], 0, psds_2.shape[0]])
plt.colorbar(label='Normalized Power')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Row in Dataset')
plt.title('Heat Map of Power Spectral Densities for Dataset 2')
plt.savefig('HP_NoP.png')


