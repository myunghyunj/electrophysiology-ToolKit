%% Batch EEG Analysis & Export
% Based on FFT_psd_histogram.m
% This script:
% 1. Processes multiple .mat files 
% 2. Saves relative power percentages to Excel files
% 3. Saves figures without annotations as .fig files

clear; close all; clc;

%% Parameters
fs = 1000; % Sampling frequency (Hz)
epoch = 2; % Epoch length in seconds

%% Select multiple .mat files
[fileNames, path] = uigetfile('*.mat', 'Select EEG data file(s)', 'MultiSelect', 'on');

% Check if user canceled file selection
if isequal(fileNames, 0)
    error('No files selected. Exiting.');
end

% Convert to cell array if only one file selected
if ~iscell(fileNames)
    fileNames = {fileNames};
end

fprintf('Processing %d files...\n', length(fileNames));

% Create output folders if they don't exist
excelFolder = fullfile(path, 'excel_results');
figureFolder = fullfile(path, 'figure_results');

if ~exist(excelFolder, 'dir')
    mkdir(excelFolder);
    fprintf('Created folder: %s\n', excelFolder);
end

if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
    fprintf('Created folder: %s\n', figureFolder);
end

%% Process each file
for fileIdx = 1:length(fileNames)
    currentFile = fileNames{fileIdx};
    fprintf('\n[%d/%d] Processing: %s\n', fileIdx, length(fileNames), currentFile);
    
    %% Load EEG Data
    try
        loadedData = load(fullfile(path, currentFile));
    catch
        fprintf('Error loading file: %s. Skipping.\n', currentFile);
        continue;
    end
    
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
    
    % If not found in 'amplifier_data', use the dynamic search
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
        fprintf('No suitable EEG data found in file: %s. Skipping.\n', currentFile);
        continue;
    end
    
    % Ensure raw_eeg_data is a row vector
    if size(raw_eeg_data,1) > 1 && size(raw_eeg_data,2) == 1
        raw_eeg_data = raw_eeg_data'; 
    end
    
    %% Create figure without annotations
    fig = figure('Position', [100, 100, 900, 700]);
    
    % Add a super title (main title) to the figure, including the filename
    % Escape underscores in the filename for proper display in the title
    title_str = sprintf('EEG Analysis: %s', strrep(currentFile, '_', '\\_'));
    sgtitle(title_str);
    
    %% EEG Bands Definition
    bands = {'Delta', [0.5,4];
             'Theta', [4,7];      % 4-7 Hz for precise theta definition
             'Alpha', [8,12];
             'Beta',  [12,30];
             'Gamma', [30,100]};
         
    num_bands = size(bands, 1);
    power_band = zeros(1,num_bands);
    
    %% Plot Raw EEG Signal
    subplot(3,1,1); % First subplot for Raw EEG
    time = (0:length(raw_eeg_data)-1) / fs;
    
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
    
    %% Analysis using bandpower function
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
    else
        gamma_theta_ratio = gamma_mean / theta_mean;
        gamma_delta_ratio = gamma_mean / delta_mean;
    end
    
    %% Original FFT Analysis
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
    
    % Calculate FFT method ratios
    if power_band(theta_idx) == 0 || power_band(delta_idx) == 0
        % Avoid division by zero
        gamma_theta_ratio_fft = NaN;
        gamma_delta_ratio_fft = NaN;
    else
        gamma_theta_ratio_fft = power_band(gamma_idx) / power_band(theta_idx);
        gamma_delta_ratio_fft = power_band(gamma_idx) / power_band(delta_idx);
    end
    
    %% Plot Frequency Histogram (without red text annotations)
    subplot(3,1,2); % Second subplot for Histogram
    bar(categorical(bands(:,1)), power_band, 'FaceColor', [0.2 0.6 0.8]);
    xlabel('EEG Frequency Bands');
    ylabel('Average Power');
    title('EEG Band Power Histogram');
    grid on;
    
    %% Plot Spectrogram
    subplot(3,1,3); % Third subplot for Spectrogram
    window_length = 300;
    overlap_length = 80;
    fft_points = 300;
    
    if length(raw_eeg_data) >= window_length
        spectrogram(raw_eeg_data, window_length, overlap_length, fft_points, fs, 'yaxis');
        ylim([0 50]);
        title('EEG Spectrogram (0-50Hz)');
    else
        text(0.5, 0.5, 'Data too short for spectrogram', 'HorizontalAlignment', 'center');
        axis off;
        title('EEG Spectrogram (Data too short)');
    end
    
    %% Save figure to .fig file
    % Remove file extension and replace with .fig
    [~, baseFileName, ~] = fileparts(currentFile);
    figFileName = fullfile(figureFolder, [baseFileName '_analysis.fig']);
    savefig(fig, figFileName);
    fprintf('Saved figure to: %s\n', figFileName);
    
    %% Save results to Excel
    % Prepare data for Excel export
    bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
    
    % Bandpower method results
    bandpowerValues = [delta_mean, theta_mean, alpha_mean, beta_mean, gamma_mean];
    bandpowerPercentages = [delta_power_rel_bp, theta_power_rel_bp, alpha_power_rel_bp, beta_power_rel_bp, gamma_power_rel_bp];
    
    % FFT method results
    fftValues = [power_band(delta_idx), power_band(theta_idx), power_band(alpha_idx), power_band(beta_idx), power_band(gamma_idx)];
    fftPercentages = [delta_power_rel, theta_power_rel, alpha_power_rel, beta_power_rel, gamma_power_rel];
    
    % Ratios
    ratios_bandpower = {'Gamma/Theta', 'Gamma/Delta'; gamma_theta_ratio, gamma_delta_ratio};
    ratios_fft = {'Gamma/Theta', 'Gamma/Delta'; gamma_theta_ratio_fft, gamma_delta_ratio_fft};
    
    % Create table for Excel
    T_bandpower = table(bandNames', bandpowerValues', bandpowerPercentages', 'VariableNames', {'Band', 'Power', 'Percentage'});
    T_fft = table(bandNames', fftValues', fftPercentages', 'VariableNames', {'Band', 'Power', 'Percentage'});
    
    % Save to Excel
    excelFileName = fullfile(excelFolder, [baseFileName '_results.xlsx']);
    
    % Write bandpower method results
    writetable(T_bandpower, excelFileName, 'Sheet', 'Bandpower Method', 'Range', 'A1');
    % Add ratios for bandpower method
    writecell({'Ratios'}, excelFileName, 'Sheet', 'Bandpower Method', 'Range', 'A8');
    writecell(ratios_bandpower, excelFileName, 'Sheet', 'Bandpower Method', 'Range', 'A9');
    
    % Write FFT method results
    writetable(T_fft, excelFileName, 'Sheet', 'FFT Method', 'Range', 'A1');
    % Add ratios for FFT method
    writecell({'Ratios'}, excelFileName, 'Sheet', 'FFT Method', 'Range', 'A8');
    writecell(ratios_fft, excelFileName, 'Sheet', 'FFT Method', 'Range', 'A9');
    
    fprintf('Saved results to: %s\n', excelFileName);
    
    % Close the figure to free memory
    close(fig);
end

fprintf('\nBatch processing complete. Processed %d files.\n', length(fileNames)); 