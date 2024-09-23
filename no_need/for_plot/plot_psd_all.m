function plot_psd_comparison()


    % Sampling information
    step_size_ms = 0.1; % step size in milliseconds
    Fs = 1000 / step_size_ms; % Sampling frequency in Hz


    csv_file_path_1 = '/Users/dracoxu/Research/data/all_raw_activities_differentp.csv';
    csv_file_path_2 = '/Users/dracoxu/Research/data/all_raw_activities_samep.csv';
    csv_file_path_3 = '/Users/dracoxu/Research/data/all_raw_activities_nop.csv';
    %csv_file_path_4 = '/Users/dracoxu/Research/data/all_raw_activities_smallp.csv';
    csv_file_path_4 = '/Users/dracoxu/Research/data/all_kkk.csv';
    

    % Load data
    data1 = load_and_process_csv(csv_file_path_1);
    data2 = load_and_process_csv(csv_file_path_2);
    data3 = load_and_process_csv(csv_file_path_3);
    data4 = load_and_process_csv(csv_file_path_4);

    % Denoise and Detrend Data
    data1 = preprocess_data(data1);
    data2 = preprocess_data(data2);
    data3 = preprocess_data(data3);
    data4 = preprocess_data(data4);

    k = 2300;
    data1 = data1(:, k+200:end);
    data2 = data2(:, k+200:end);
    data3 = data3(:, k+200:end);
    data4 = data4(:, k+200:end);
    




    % Intermediate Steps Visualization (After Preprocessing)
    figure;
    subplot(4, 1, 1);
    plot(mean(data1, 1), 'b', 'LineWidth', 2);
    title('Example Signal After Preprocessing for one trial in Dataset 1');
    xlabel('Sample Number');
    ylabel('Amplitude');

    subplot(4, 1, 2);
    plot(mean(data2, 1), 'r', 'LineWidth', 2);
    title('Example Signal After Preprocessing for one trial in Dataset 2');
    xlabel('Sample Number');
    ylabel('Amplitude');

    subplot(4, 1, 3);
    plot(mean(data3, 1), 'k', 'LineWidth', 2);
    title('Example Signal After Preprocessing for one trial in Dataset 3');
    xlabel('Sample Number');
    ylabel('Amplitude');

    subplot(4, 1, 4);
    plot(mean(data4, 1), 'g', 'LineWidth', 2);
    title('Example Signal After Preprocessing for one trial in Dataset 4');
    xlabel('Sample Number');
    ylabel('Amplitude');

    % Compute average PSD
    [psd1, freq1] = average_psd(data1, Fs);
    [psd2, freq2] = average_psd(data2, Fs);
    [psd3, freq3] = average_psd(data3, Fs);
    [psd4, freq4] = average_psd(data4, Fs);
    



        % Plot the results of PSDs
    figure;
    plot(freq1(1:25), psd1(1:25), 'b', 'LineWidth', 2);
    hold on;
    plot(freq2(1:25), psd2(1:25), 'r', 'LineWidth', 2);
    hold on;
    plot(freq3(1:25), psd3(1:25), 'k', 'LineWidth', 2);
    hold on;
    plot(freq4(1:25), psd4(1:25), 'g', 'LineWidth', 2);
    title('Unnormalized Average Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power');
    legend('Diff P', 'Same P','No P','Small P');
    grid on;
    hold off;


    z = 3
    psd1= psd1(z:end);
    psd2= psd2(z:end);
    psd3= psd3(z:end);
    psd4= psd4(z:end);
    freq1 =freq1(z:end);
    freq2 =freq2(z:end);
    freq3 =freq3(z:end);
    freq4 =freq4(z:end);

    % Plot the results of PSDs with spline interpolation for smoothing
figure;

f_min = 3.40;  % Set the minimum frequency for interpolation
f_max = max(freq1(1:25));  % Maximum frequency based on your data
f_interp = linspace(f_min, f_max, 100);  % Create linear space starting from 5 Hz

% Interpolate PSDs using the new frequency range starting from 5 Hz
psd1_interp = spline(freq1(1:25), psd1(1:25), f_interp);  % Spline interpolation for dataset 1
psd2_interp = spline(freq2(1:25), psd2(1:25), f_interp);  % Spline interpolation for dataset 2
psd3_interp = spline(freq3(1:25), psd3(1:25), f_interp);  % Spline interpolation for dataset 2
psd4_interp = spline(freq4(1:25), psd4(1:25), f_interp);  % Spline interpolation for dataset 2

    psd1_interp = psd1_interp / max(psd1_interp);
    psd2_interp = psd2_interp / max(psd2_interp);
    psd3_interp = psd3_interp / max(psd3_interp);
    psd4_interp = psd4_interp / max(psd4_interp);




plot(f_interp, psd1_interp, 'b', 'LineWidth', 2);
hold on;
plot(f_interp, psd2_interp, 'r', 'LineWidth', 2);
hold on;
plot(f_interp, psd3_interp, 'k', 'LineWidth', 2);
hold on;
plot(f_interp, psd4_interp, 'g', 'LineWidth', 2);
title('normalized Average Power Spectral Density');
xlabel('Frequency (Hz)');
xlim([3 35]); % Set x-axis limits starting from 5 Hz
ylabel('Power');
legend('Diff P', 'Same P','No P','Small P');
grid on;
hold off;









    % Normalize the PSDs so that their max is 1
    psd1 = psd1 / max(psd1);
    psd2 = psd2 / max(psd2);
    psd3 = psd3 / max(psd3);
    psd4 = psd4 / max(psd4);

    % Plot the results of PSDs
    figure;
    plot(freq1(1:25), psd1(1:25), 'b', 'LineWidth', 2);
    hold on;
    plot(freq2(1:25), psd2(1:25), 'r', 'LineWidth', 2);
    hold on;
    plot(freq3(1:25), psd3(1:25), 'k', 'LineWidth', 2);
    hold on;
    plot(freq4(1:25), psd4(1:25), 'g', 'LineWidth', 2);
    title('Normalized Average Power Spectral Density');
    xlabel('Frequency (Hz)');
    ylabel('Power');
    legend('Diff P', 'Same P','No P','Small P');
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
    window_size = 0.65 * length(data(1, :));  % Using the full data length of the first trial
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