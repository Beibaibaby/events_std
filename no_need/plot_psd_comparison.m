function plot_psd_comparison()


    % Sampling information
    step_size_ms = 0.1; % step size in milliseconds
    Fs = 1000 / step_size_ms; % Sampling frequency in Hz


    csv_file_path_1 = '/Users/dracoxu/Research/data/all_raw_activities_differentp.csv';
    csv_file_path_2 = '/Users/dracoxu/Research/data/all_raw_activities_samep.csv';

    % Load data
    data1 = load_and_process_csv(csv_file_path_1);
    data2 = load_and_process_csv(csv_file_path_2);

    % Denoise and Detrend Data
    data1 = preprocess_data(data1);
    data2 = preprocess_data(data2);

    k = 2300;
    data1 = data1(:, k+200:end);
    data2 = data2(:, k+200:end);




    % Intermediate Steps Visualization (After Preprocessing)
    figure;
    subplot(2, 1, 1);
    plot(mean(data1, 1), 'b', 'LineWidth', 2);
    title('Example Signal After Preprocessing for one trial in Dataset 1');
    xlabel('Sample Number');
    ylabel('Amplitude');

    subplot(2, 1, 2);
    plot(mean(data2, 1), 'r', 'LineWidth', 2);
    title('Example Signal After Preprocessing for one trial in Dataset 2');
    xlabel('Sample Number');
    ylabel('Amplitude');

    % Compute average PSD
    [psd1, freq1] = average_psd(data1, Fs);
    [psd2, freq2] = average_psd(data2, Fs);
    



        % Plot the results of PSDs
    figure;
    plot(freq1(1:25), psd1(1:25), 'b', 'LineWidth', 2);
    hold on;
    plot(freq2(1:25), psd2(1:25), 'r', 'LineWidth', 2);
    title('Unnormalized Average Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power');
    legend('Diff P', 'Same P');
    grid on;
    hold off;


    z = 3
    psd1= psd1(z:end);
    psd2= psd2(z:end);
    freq1 =freq1(z:end);
    freq2 =freq2(z:end);

    % Plot the results of PSDs with spline interpolation for smoothing
figure;
f_interp = linspace(min(freq1(1:25)), max(freq1(1:25)), 500);  % increase the number of points in the frequency range for a smooth plot
psd1_interp = spline(freq1(1:25), psd1(1:25), f_interp);  % spline interpolation
psd2_interp = spline(freq2(1:25), psd2(1:25), f_interp);  % spline interpolation


f_min = 5;  % Set the minimum frequency for interpolation
f_max = max(freq1(1:25));  % Maximum frequency based on your data
f_interp = linspace(f_min, f_max, 1000);  % Create linear space starting from 5 Hz

% Interpolate PSDs using the new frequency range starting from 5 Hz
psd1_interp = spline(freq1(1:25), psd1(1:25), f_interp);  % Spline interpolation for dataset 1
psd2_interp = spline(freq2(1:25), psd2(1:25), f_interp);  % Spline interpolation for dataset 2

    psd1_interp = psd1_interp / max(psd1_interp);
    psd2_interp = psd2_interp / max(psd2_interp);




plot(f_interp, psd1_interp, 'b', 'LineWidth', 2);
hold on;
plot(f_interp, psd2_interp, 'r', 'LineWidth', 2);
title('normalized Average Power Spectral Density');
xlabel('Frequency (Hz)');
xlim([5 35]); % Set x-axis limits starting from 5 Hz
ylabel('Power');
legend('Diff P', 'Same P');
grid on;
hold off;









    % Normalize the PSDs so that their max is 1
    psd1 = psd1 / max(psd1);
    psd2 = psd2 / max(psd2);

    % Plot the results of PSDs
    figure;
    plot(freq1(1:25), psd1(1:25), 'b', 'LineWidth', 2);
    hold on;
    plot(freq2(1:25), psd2(1:25), 'r', 'LineWidth', 2);
    title('Normalized Average Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power');
    legend('Diff P', 'Same P');
    grid on;
    hold off;

end

function matrix = load_and_process_csv(csv_file_path)
    % Load CSV file
    matrix = readmatrix(csv_file_path);
end

function data = preprocess_data(data)
    % Apply Savitzky-Golay filter to denoise data
    for i = 1:size(data, 1)
        data(i, :) = sgolayfilt(data(i, :), 1, 55);  % first-order polynomial, window size of 49
    end
    
     %Center data by subtracting the mean of each trial
    %data = data - mean(data, 2);
    %for i = 1:size(data,1)
    %    data(i,:) = detrend(data(1,:));
    %end

end

function data = preprocess_data_with_detrend(data)
    % Apply Savitzky-Golay filter to denoise data
    for i = 1:size(data, 1)
        data(i, :) = sgolayfilt(data(i, :), 1, 209);  % first-order polynomial, window size of 49
    end
    
     
    
    for i = 1:size(data,1)
        data(i,:) = detrend(data(1,:));
    end

end

function [avg_psd, freq] = average_psd(data, Fs)
    % Specifying a larger window size to increase frequency resolution
    window_size = 0.35 * length(data(1, :));  % Using the full data length of the first trial
    window = hamming(window_size);  % Using a Hamming window to reduce spectral leakage

    [~, f] = pwelch(data(1, :), window, [], [], Fs);
    psds = zeros(size(data, 1), length(f));

    for i = 1:size(data, 1)
        [psd, ~] = pwelch(data(i, :), window, [], [], Fs);
        psds(i, :) = psd.';  % Transpose psd to match dimensions
    end

    avg_psd = mean(psds, 1);
    freq = f;
end


function matrix = load_and_expand_csv(csv_file_path)
    % Load CSV file
    original_matrix = readmatrix(csv_file_path);
    % Expand each column 10 times
    matrix = repelem(original_matrix, 1, 10);
end