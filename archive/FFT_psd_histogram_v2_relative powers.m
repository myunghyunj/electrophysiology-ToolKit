%% EEG Analysis & Frequency Histogram Plotting (Robust Data Loader)

clear; close all; clc;

%% Parameters
fs = 1000; % Sampling frequency (Hz)
epoch = 2; % Epoch length in seconds

%% Load EEG Data (Robust)
[file, path] = uigetfile('*.mat', 'Select EEG data file');
loadedData = load(fullfile(path, file));

raw_eeg_data = []; % Initialize raw_eeg_data
found_data = false;

% Attempt to load from 'amplifier_data' field first
if isfield(loadedData, 'amplifier_data')
    potential_field_data = loadedData.amplifier_data;
    if isnumeric(potential_field_data) && size(potential_field_data, 1) >= 1
        % Try first row of amplifier_data
        temp_eeg = potential_field_data(1,:);
        if length(temp_eeg) > 1 && sum(abs(temp_eeg)) > 1e-9 % Check for non-zero/non-empty and non-scalar signal
            raw_eeg_data = temp_eeg;
            fprintf('EEG data loaded from field: amplifier_data, row 1\n');
            found_data = true;
        elseif size(potential_field_data, 1) >= 2
            % First row was unsuitable, try second row of amplifier_data
            temp_eeg = potential_field_data(2,:);
            if length(temp_eeg) > 1 && sum(abs(temp_eeg)) > 1e-9 % Check for non-zero/non-empty and non-scalar signal
                raw_eeg_data = temp_eeg;
                fprintf('EEG data loaded from field: amplifier_data, row 2 (row 1 was empty/zero/scalar)\n');
                found_data = true;
            end
        end
    end
end

% If not found in 'amplifier_data', use the original dynamic search
if ~found_data
    fields = fieldnames(loadedData);
    for i = 1:length(fields)
        current_field_data = loadedData.(fields{i});
        if isnumeric(current_field_data) && size(current_field_data,1) >= 1
            % Try first row
            temp_eeg = current_field_data(1,:);
            if length(temp_eeg) > 1 && sum(abs(temp_eeg)) > 1e-9 % Check for non-zero/non-empty and non-scalar signal
                raw_eeg_data = temp_eeg;
                fprintf('EEG data loaded from field: %s, row 1\n', fields{i});
                found_data = true;
                break;
            elseif size(current_field_data,1) >= 2
                % First row was unsuitable, try second row
                temp_eeg = current_field_data(2,:);
                if length(temp_eeg) > 1 && sum(abs(temp_eeg)) > 1e-9 % Check for non-zero/non-empty and non-scalar signal
                    raw_eeg_data = temp_eeg;
                    fprintf('EEG data loaded from field: %s, row 2 (row 1 was empty/zero/scalar)\n', fields{i});
                    found_data = true;
                    break;
                end
            end
        end
    end
end

if ~found_data || isempty(raw_eeg_data)
    error('No suitable EEG data (non-empty row in a numeric field) found in the loaded .mat file.');
end

% Ensure raw_eeg_data is a row vector
if size(raw_eeg_data,1) > 1 && size(raw_eeg_data,2) == 1
    raw_eeg_data = raw_eeg_data'; 
end

%% Raw EEG Signal Plot
time = (0:length(raw_eeg_data)-1) / fs;
figure; % This creates the single figure window for all subplots

% Add a super title (main title) to the figure, including the filename
% Escape underscores in the filename for proper display in the title
if exist('file', 'var') % Check if file variable exists (it should from uigetfile)
    title_str = sprintf('EEG Analysis: %s', strrep(file, '_', '\\_'));
else
    title_str = 'EEG Analysis';
end
sgtitle(title_str);

subplot(3,1,1); % First subplot for Raw EEG
if length(time) < 2
    plot(time, raw_eeg_data, 'o-'); % Mark points if few
    title('Raw EEG Signal (Short Data)');
    if isempty(time)
        xlim_vals = [0 1/fs]; % Default if time is totally empty
    elseif min(time) == max(time)
        xlim_vals = [time(1)-0.5/fs, time(1)+0.5/fs]; % Give some space for single point
    else
        % Should not happen if length(time) < 2 and not empty and not single point
        xlim_vals = [min(time), max(time) + eps]; % Fallback, ensure range
    end
    warning('EEG data is too short for a standard plot. Displaying available points.');
else
    plot(time, raw_eeg_data);
    title('Raw EEG Signal');
    xlim_vals = [0, max(time)];
end

xlabel('Time (s)');
ylabel('EEG Amplitude');
grid on;
if ~isempty(xlim_vals) % Ensure xlim_vals is set
    xlim(xlim_vals);
else
    % Fallback xlim if something went wrong with xlim_vals logic
    if ~isempty(time)
        xlim([min(time), max(time)+eps]);
    else
        xlim([0,1]); % Absolute fallback
    end
end

%% EEG Bands Definition
bands = {'Delta', [0.5,4];
         'Theta', [4,7];      % 4-7 Hz for precise theta definition
         'Alpha', [8,12];
         'Beta',  [12,30];
         'Gamma', [30,100]};
     
num_bands = size(bands, 1);
power_band = zeros(1,num_bands);

%% Hypothesis Test: Enhanced analysis using bandpower function
data_duration = length(raw_eeg_data)/fs;
epoch_length = 2; % 2-second epochs, ideal for delta/theta analysis
num_epochs_bandpower = floor(data_duration/epoch_length);
epoch_samples_bandpower = epoch_length*fs;

% Initialize power arrays for each band
gamma_power = zeros(1, num_epochs_bandpower);
theta_power = zeros(1, num_epochs_bandpower);
delta_power = zeros(1, num_epochs_bandpower);
alpha_power = zeros(1, num_epochs_bandpower);
beta_power = zeros(1, num_epochs_bandpower);

for ep = 1:num_epochs_bandpower
    % Extract current epoch
    start_idx = (ep-1)*epoch_samples_bandpower + 1;
    end_idx = min(ep*epoch_samples_bandpower, length(raw_eeg_data)); % Ensure we don't exceed data length
    if end_idx - start_idx + 1 < epoch_samples_bandpower/2
        % Skip if epoch is less than half of desired length
        continue;
    end
    epoch_data = raw_eeg_data(start_idx:end_idx);
    
    % Calculate band power using bandpower function (more accurate)
    try
        % Gamma (30-100 Hz)
        gamma_power(ep) = bandpower(epoch_data, fs, [30 100]);
        % Beta (12-30 Hz)
        beta_power(ep) = bandpower(epoch_data, fs, [12 30]);
        % Alpha (8-12 Hz)
        alpha_power(ep) = bandpower(epoch_data, fs, [8 12]);
        % Theta (4-7 Hz)
        theta_power(ep) = bandpower(epoch_data, fs, [4 7]);
        % Delta (0.5-4 Hz)
        delta_power(ep) = bandpower(epoch_data, fs, [0.5 4]);
    catch
        % Fallback if bandpower is not available (older MATLAB versions)
        fprintf('Warning: bandpower function not available, using manual calculation for epoch %d\n', ep);
        % Manual calculation of band power using FFT
        L = length(epoch_data);
        Y = fft(epoch_data);
        P2 = abs(Y/L);
        P1 = P2(1:floor(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        freq = fs*(0:(floor(L/2)))/L;
        
        % Manually find indices for each band
        gamma_idx = (freq >= 30) & (freq <= 100);
        beta_idx = (freq >= 12) & (freq <= 30);
        alpha_idx = (freq >= 8) & (freq <= 12);
        theta_idx = (freq >= 4) & (freq <= 7);
        delta_idx = (freq >= 0.5) & (freq <= 4);
        
        gamma_power(ep) = sum(P1(gamma_idx).^2);
        beta_power(ep) = sum(P1(beta_idx).^2);
        alpha_power(ep) = sum(P1(alpha_idx).^2);
        theta_power(ep) = sum(P1(theta_idx).^2);
        delta_power(ep) = sum(P1(delta_idx).^2);
    end
end

% Calculate mean powers (removing any zeros from skipped epochs)
gamma_mean = mean(gamma_power(gamma_power > 0));
beta_mean = mean(beta_power(beta_power > 0));
alpha_mean = mean(alpha_power(alpha_power > 0));
theta_mean = mean(theta_power(theta_power > 0));
delta_mean = mean(delta_power(delta_power > 0));

if isempty(gamma_mean), gamma_mean = 0; end
if isempty(beta_mean), beta_mean = 0; end
if isempty(alpha_mean), alpha_mean = 0; end
if isempty(theta_mean), theta_mean = 0; end
if isempty(delta_mean), delta_mean = 0; end

% Calculate total power for relative power calculations
total_power_bandpower = gamma_mean + beta_mean + alpha_mean + theta_mean + delta_mean;

% Calculate relative powers (percentages)
if total_power_bandpower > 0
    delta_power_rel_bp = (delta_mean / total_power_bandpower) * 100;
    theta_power_rel_bp = (theta_mean / total_power_bandpower) * 100;
    alpha_power_rel_bp = (alpha_mean / total_power_bandpower) * 100;
    beta_power_rel_bp = (beta_mean / total_power_bandpower) * 100;
    gamma_power_rel_bp = (gamma_mean / total_power_bandpower) * 100;
else
    delta_power_rel_bp = 0;
    theta_power_rel_bp = 0;
    alpha_power_rel_bp = 0;
    beta_power_rel_bp = 0;
    gamma_power_rel_bp = 0;
end

% Calculate ratios with error handling
if theta_mean == 0 || delta_mean == 0
    % Avoid division by zero
    gamma_theta_ratio = NaN;
    gamma_delta_ratio = NaN;
    ratio_str_theta = 'N/A (Theta power is 0)';
    ratio_str_delta = 'N/A (Delta power is 0)';
else
    gamma_theta_ratio = gamma_mean / theta_mean;
    gamma_delta_ratio = gamma_mean / delta_mean;
    ratio_str_theta = sprintf('%.2f', gamma_theta_ratio);
    ratio_str_delta = sprintf('%.2f', gamma_delta_ratio);
end

% Display results in command window
fprintf('Hypothesis Test Results (Bandpower Method):\n');
fprintf('Gamma Power: %.4e (%.2f%%)\n', gamma_mean, gamma_power_rel_bp);
fprintf('Beta Power: %.4e (%.2f%%)\n', beta_mean, beta_power_rel_bp);
fprintf('Alpha Power: %.4e (%.2f%%)\n', alpha_mean, alpha_power_rel_bp);
fprintf('Theta Power: %.4e (%.2f%%)\n', theta_mean, theta_power_rel_bp);
fprintf('Delta Power: %.4e (%.2f%%)\n', delta_mean, delta_power_rel_bp);
fprintf('Gamma/Theta Ratio: %s\n', ratio_str_theta);
fprintf('Gamma/Delta Ratio: %s\n', ratio_str_delta);

% Store results for later use if needed
result_data.gamma_power = gamma_mean;
result_data.beta_power = beta_mean;
result_data.alpha_power = alpha_mean;
result_data.theta_power = theta_mean;
result_data.delta_power = delta_mean;
result_data.gamma_theta_ratio = gamma_theta_ratio;
result_data.gamma_delta_ratio = gamma_delta_ratio;
result_data.delta_power_rel = delta_power_rel_bp;
result_data.theta_power_rel = theta_power_rel_bp;
result_data.alpha_power_rel = alpha_power_rel_bp;

%% Original FFT Analysis (we keep this for backward compatibility)
data_duration = length(raw_eeg_data)/fs;
num_epochs = floor(data_duration/epoch);
epoch_samples = epoch*fs;

for ep = 1:num_epochs
    epoch_data = raw_eeg_data((ep-1)*epoch_samples+1 : ep*epoch_samples);
    
    % FFT computation
    L = length(epoch_data);
    Y = fft(epoch_data);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    freq = fs*(0:(L/2))/L;

    % Accumulate power per band
    for b = 1:num_bands
        freq_range = bands{b,2};
        band_indices = freq >= freq_range(1) & freq <= freq_range(2);
        power_band(b) = power_band(b) + sum(P1(band_indices).^2);
    end
end

% Average power across epochs
power_band = power_band / num_epochs;

% Calculate total power for relative power calculations
total_power = sum(power_band);

% Calculate relative powers (percentages)
delta_idx = find(strcmp(bands(:,1), 'Delta'));
theta_idx = find(strcmp(bands(:,1), 'Theta'));
alpha_idx = find(strcmp(bands(:,1), 'Alpha'));
beta_idx = find(strcmp(bands(:,1), 'Beta'));
gamma_idx = find(strcmp(bands(:,1), 'Gamma'));

if total_power > 0
    delta_power_rel = (power_band(delta_idx) / total_power) * 100;
    theta_power_rel = (power_band(theta_idx) / total_power) * 100;
    alpha_power_rel = (power_band(alpha_idx) / total_power) * 100;
    beta_power_rel = (power_band(beta_idx) / total_power) * 100;
    gamma_power_rel = (power_band(gamma_idx) / total_power) * 100;
else
    delta_power_rel = 0;
    theta_power_rel = 0;
    alpha_power_rel = 0;
    beta_power_rel = 0;
    gamma_power_rel = 0;
end

% Display relative powers
fprintf('\nRelative Powers (Original FFT Method):\n');
fprintf('Delta power: %.2f%%\n', delta_power_rel);
fprintf('Theta power: %.2f%%\n', theta_power_rel);
fprintf('Alpha power: %.2f%%\n', alpha_power_rel);
fprintf('Beta power: %.2f%%\n', beta_power_rel);
fprintf('Gamma power: %.2f%%\n', gamma_power_rel);

%% Plot Frequency Histogram
subplot(3,1,2); % Second subplot for Histogram
bar(categorical(bands(:,1)), power_band, 'FaceColor', [0.2 0.6 0.8]);
xlabel('EEG Frequency Bands');
ylabel('Average Power');
title('EEG Band Power Histogram');
grid on;

% Calculate Gamma, Delta power and their ratio for display
if ~isempty(gamma_idx) && ~isempty(delta_idx) && ~isempty(theta_idx) && length(power_band) >= max([gamma_idx, delta_idx, theta_idx])
    gamma_power_orig = power_band(gamma_idx);
    delta_power_orig = power_band(delta_idx);
    theta_power_orig = power_band(theta_idx);
    alpha_power_orig = power_band(alpha_idx);

    % For the original FFT method values
    if delta_power_orig == 0
        gamma_delta_ratio_str = 'N/A (Delta power is 0)';
    else
        gamma_delta_ratio_orig = gamma_power_orig / delta_power_orig;
        gamma_delta_ratio_str = sprintf('%.2f', gamma_delta_ratio_orig);
    end
    
    % Add calculation for Gamma/Theta ratio from original FFT
    if theta_power_orig == 0
        gamma_theta_ratio_str = 'N/A (Theta power is 0)';
    else
        gamma_theta_ratio_orig = gamma_power_orig / theta_power_orig;
        gamma_theta_ratio_str = sprintf('%.2f', gamma_theta_ratio_orig);
    end

    % Prepare display text with both bandpower and original FFT results
    text_to_display = {
        '--- Bandpower Method ---', ...
        sprintf('Gamma Power: %.2e (%.1f%%)', gamma_mean, gamma_power_rel_bp), ...
        sprintf('Beta Power: %.2e (%.1f%%)', beta_mean, beta_power_rel_bp), ...
        sprintf('Alpha Power: %.2e (%.1f%%)', alpha_mean, alpha_power_rel_bp), ...
        sprintf('Theta Power: %.2e (%.1f%%)', theta_mean, theta_power_rel_bp), ...
        sprintf('Delta Power: %.2e (%.1f%%)', delta_mean, delta_power_rel_bp), ...
        sprintf('Gamma/Theta Ratio: %s', ratio_str_theta), ...
        sprintf('Gamma/Delta Ratio: %s', ratio_str_delta), ...
        '', ...
        '--- Original FFT Method ---', ...
        sprintf('Gamma Power: %.2f (%.1f%%)', gamma_power_orig, gamma_power_rel), ...
        sprintf('Beta Power: %.2f (%.1f%%)', power_band(beta_idx), beta_power_rel), ...
        sprintf('Alpha Power: %.2f (%.1f%%)', alpha_power_orig, alpha_power_rel), ...
        sprintf('Theta Power: %.2f (%.1f%%)', theta_power_orig, theta_power_rel), ...
        sprintf('Delta Power: %.2f (%.1f%%)', delta_power_orig, delta_power_rel), ...
        sprintf('Gamma/Theta Ratio: %s', gamma_theta_ratio_str), ...
        sprintf('Gamma/Delta Ratio: %s', gamma_delta_ratio_str)
    };

    % Get current axis to position text
    ax = gca;
    % Ensure XLim and YLim are treated as numeric for range calculation
    current_xlim_numeric = double(ax.XLim);
    current_ylim_numeric = double(ax.YLim);

    % We position slightly to the right of the very edge of the plot area.
    text_x_position = current_xlim_numeric(1) + 0.1 * (current_xlim_numeric(2) - current_xlim_numeric(1)); 
    % Position text near the top of the y-axis.
    text_y_position = current_ylim_numeric(2) - 0.05 * (current_ylim_numeric(2) - current_ylim_numeric(1)); 

    text(text_x_position, text_y_position, text_to_display, ...
        'Color', 'red', ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left');
else
    disp('Could not find Gamma, Delta, or Theta bands, or power_band array is too short to display ratio.');
end

%% Spectrogram (Optional Visualization)
subplot(3,1,3); % Third subplot for Spectrogram
window_length = 300;
overlap_length = 80;
fft_points = 300;

if length(raw_eeg_data) >= window_length
    spectrogram(raw_eeg_data, window_length, overlap_length, fft_points, fs, 'yaxis');
    ylim([0 50]);
    title('EEG Spectrogram (0-50Hz)');
else
    disp('Raw EEG data is too short to generate a spectrogram with the current settings.');
    title('EEG Spectrogram (Data too short)');
    % Optionally, you could try to plot something simpler or nothing at all
    % For example, just show an empty axes or a text message
    text(0.5, 0.5, 'Data too short for spectrogram', 'HorizontalAlignment', 'center');
    axis off;
end 