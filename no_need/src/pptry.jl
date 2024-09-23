ENV["GKSwstype"] = "100"  # Use the off-screen GKS terminal
using Plots
using DSP  # For signal processing functions
using FFTW  # For FFT functions

# Your data
data = [1.39, 1.4400000000000002, 1.35, 1.4, 1.4, 1.5, 2.12, 3.76, 5.58, 6.43, 9.520000000000001, 10.040000000000001, 11.620000000000001, 9.37, 7.1000000000000005, 5.0600000000000005, 2.19, 0.8, 0.95, 1.86, 1.85, 1.97, 2.39, 2.29, 1.79, 2.78, 2.97, 2.84, 3.3899999999999997, 3.14, 2.83, 3.0599999999999996, 3.18, 3.0799999999999996, 3.16, 3.12, 2.82, 2.32, 1.9200000000000002, 2.06, 2.18, 2.52, 2.73, 2.75, 2.97, 2.7399999999999998, 2.61, 2.64, 2.27, 1.58, 1.6199999999999999, 1.6, 2.07, 2.0799999999999996, 2.44, 2.19, 2.26, 1.58, 1.74, 1.2899999999999998, 1.1900000000000002, 1.52, 1.45, 1.22, 1.55, 1.9, 2.09, 2.33, 2.15, 2.1, 1.9, 1.5399999999999998, 1.9, 1.82, 1.75, 1.64, 1.61, 1.61, 1.8699999999999999, 1.8, 1.6900000000000002, 1.57, 1.32, 1.17, 1.11, 1.36, 1.1800000000000002, 1.27, 1.7, 1.74, 1.65, 1.74, 1.52, 1.39, 1.5299999999999998, 1.55, 1.79, 2.0799999999999996, 1.74, 1.73, 1.9400000000000002]

fs = 1 / 0.020  # 50 Hz, as each bin is 20ms
nperseg = length(data)

# Hann window for Julia. Note: In Julia, it might be called differently, ensure to use the correct window function.
hann_window = hanning(nperseg)

# Compute the power spectral density using Welch's method
psd = periodogram(data, fs=fs, window=hann_window, nfft=nperseg)

# Normalize the power spectral density
normalized_psd = psd.power / sum(psd.power)

# Plotting
p=plot(psd.freq, normalized_psd, label="Normalized PSD", xlabel="Frequency (Hz)", ylabel="Normalized Power", title="Normalized Power Spectral Density")

savefig(p, "normalized_psd.png")