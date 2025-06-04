%% EEG Analysis & Frequency Histogram Plotting (Robust Data Loader)

clear; close all; clc;

%% Parameters
fs = 1000; % Sampling frequency (Hz)
epoch = 1; % Epoch length in seconds

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
         'Theta', [4,8];
         'Alpha', [8,12];
         'Beta',  [12,30];
         'Gamma', [30,100]};
     
num_bands = size(bands, 1);
power_band = zeros(1,num_bands);

%% Epoch and FFT Analysis
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

%% Plot Frequency Histogram
subplot(3,1,2); % Second subplot for Histogram
bar(categorical(bands(:,1)), power_band, 'FaceColor', [0.2 0.6 0.8]);
xlabel('EEG Frequency Bands');
ylabel('Average Power');
title('EEG Band Power Histogram');
grid on;

% Calculate Gamma, Delta power and their ratio
gamma_idx = find(strcmp(bands(:,1), 'Gamma'));
delta_idx = find(strcmp(bands(:,1), 'Delta'));

if ~isempty(gamma_idx) && ~isempty(delta_idx) && length(power_band) >= max(gamma_idx, delta_idx)
    gamma_power = power_band(gamma_idx);
    delta_power = power_band(delta_idx);

    if delta_power == 0
        gamma_delta_ratio_str = 'N/A (Delta power is 0)';
    else
        gamma_delta_ratio = gamma_power / delta_power;
        gamma_delta_ratio_str = sprintf('%.2f', gamma_delta_ratio);
    end

    text_to_display = {
        sprintf('Gamma Power: %.2f', gamma_power), ...
        sprintf('Delta Power: %.2f', delta_power), ...
        ['Gamma/Delta Ratio: ', gamma_delta_ratio_str]
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
    disp('Could not find Gamma or Delta bands, or power_band array is too short to display ratio.');
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
