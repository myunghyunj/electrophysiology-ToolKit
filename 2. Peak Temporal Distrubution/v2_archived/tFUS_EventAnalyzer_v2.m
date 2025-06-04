%% tFUS_EventAnalyzer_v2.m
% MATLAB script for EEG event analysis and visualization
% This script:
% 1. Processes multiple .mat files containing EEG recording data
% 2. Identifies peaks using robust RMS-based spike detection
% 3. Groups detected peaks by condition (STIM vs SHAM)
% 4. Creates visualization of peak distribution and temporal dynamics

clear; close all; clc;

%% Analysis Parameters
params.win_sec          = 0.2;   % RMS window (s) for EEG event detection
params.k_thresh         = 3.5;   % Multiplier for EEG RMS threshold (baseline = median(RMS))
params.min_dur_ms       = 50;    % Minimum event duration for EEG (ms)
params.refractory_ms    = 50;    % Refractory period to merge EEG events (ms)

% Set up output directory
outputDir = fullfile(pwd, 'analysis', 'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Folder selection for STIM and SHAM data
% Ask user to select STIM folder
stim_pathname = uigetdir(pwd, 'Select folder containing STIM .mat files');
if isequal(stim_pathname, 0)
    error('No STIM folder selected. Aborting.');
end

% Get all .mat files in the STIM folder
stim_files = dir(fullfile(stim_pathname, '*.mat'));
stim_filenames = {stim_files.name};
num_stim_files = length(stim_filenames);

if num_stim_files == 0
    error('No .mat files found in the selected STIM folder.');
end
fprintf('Found %d .mat files in the STIM folder.\n', num_stim_files);

% Ask user to select SHAM folder
sham_pathname = uigetdir(pwd, 'Select folder containing SHAM .mat files');
if isequal(sham_pathname, 0)
    error('No SHAM folder selected. Aborting.');
end

% Get all .mat files in the SHAM folder
sham_files = dir(fullfile(sham_pathname, '*.mat'));
sham_filenames = {sham_files.name};
num_sham_files = length(sham_filenames);

if num_sham_files == 0
    error('No .mat files found in the selected SHAM folder.');
end
fprintf('Found %d .mat files in the SHAM folder.\n', num_sham_files);

% Get all user parameters in one dialog (except file inputs)
prompt = {'Enter bin size in minutes (e.g., 1 or 5):'};
dlgtitle = 'Analysis Parameters';
dims = [1 50];
definput = {'5'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

% Process bin size
if isempty(answer) || isempty(answer{1})
    binMin = 5; % Default to 5 minutes if cancelled
    fprintf('Using default bin size of 5 minutes.\n');
else
    binMin = str2double(answer{1});
    if isnan(binMin) || binMin <= 0
        binMin = 5; % Default to 5 minutes if invalid input
        fprintf('Invalid bin size entered. Using default bin size of 5 minutes.\n');
    else
        fprintf('Using bin size of %g minutes.\n', binMin);
    end
end

% Initialize arrays to store peak data
peak_time_STIM = [];
z_STIM = [];
filenames_STIM = {};  % Track source filename for each peak
peak_time_SHAM = [];
z_SHAM = [];
filenames_SHAM = {};  % Track source filename for each peak

%% Process STIM files
fprintf('\nProcessing STIM files...\n');

% Loop through each STIM file
for i = 1:num_stim_files
    current_file = fullfile(stim_pathname, stim_filenames{i});
    fprintf('Processing STIM file %d/%d: %s\n', i, num_stim_files, stim_filenames{i});
    
    try
        % Run actual peak detection
        [file_peaks, file_z] = detectPeaksFromFile(current_file, params);
        
        if ~isempty(file_peaks)
            % Add to collection
            peak_time_STIM = [peak_time_STIM, file_peaks];
            z_STIM = [z_STIM, file_z];
            % Track filename for each peak
            filenames_STIM = [filenames_STIM, repmat(stim_filenames(i), 1, length(file_peaks))];
            fprintf('  Found %d peaks\n', length(file_peaks));
        else
            fprintf('  No peaks detected\n');
        end
    catch err
        warning('Error processing file %s: %s', stim_filenames{i}, err.message);
    end
end

%% Process SHAM files
fprintf('\nProcessing SHAM files...\n');

% Loop through each SHAM file
for i = 1:num_sham_files
    current_file = fullfile(sham_pathname, sham_filenames{i});
    fprintf('Processing SHAM file %d/%d: %s\n', i, num_sham_files, sham_filenames{i});
    
    try
        % Run actual peak detection
        [file_peaks, file_z] = detectPeaksFromFile(current_file, params);
        
        if ~isempty(file_peaks)
            % Add to collection
            peak_time_SHAM = [peak_time_SHAM, file_peaks];
            z_SHAM = [z_SHAM, file_z];
            % Track filename for each peak
            filenames_SHAM = [filenames_SHAM, repmat(sham_filenames(i), 1, length(file_peaks))];
            fprintf('  Found %d peaks\n', length(file_peaks));
        else
            fprintf('  No peaks detected\n');
        end
    catch err
        warning('Error processing file %s: %s', sham_filenames{i}, err.message);
    end
end

%% Sort data by time 
if ~isempty(peak_time_STIM)
    [peak_time_STIM, idx] = sort(peak_time_STIM);
    z_STIM = z_STIM(idx);
    filenames_STIM = filenames_STIM(idx);
end

if ~isempty(peak_time_SHAM)
    [peak_time_SHAM, idx] = sort(peak_time_SHAM);
    z_SHAM = z_SHAM(idx);
    filenames_SHAM = filenames_SHAM(idx);
end

%% Generate temporal PSD map with user-specified bin size
fprintf('\nPlotting temporal PSD map with %g-minute bins...\n', binMin);

% Generate the temporal PSD map
plotTemporalPSDmap(peak_time_STIM, z_STIM, ...
                   peak_time_SHAM, z_SHAM, ...
                   30, binMin, outputDir, ...
                   filenames_STIM, filenames_SHAM, ...
                   num_stim_files, num_sham_files);




%% Create Figure: Z-score Distribution Bar Plot
fprintf('\nGenerating Z-score distribution bar plot...\n');
fig_zscore = figure('Position',[100 100 800 600],'Name','Z-score Distribution Bar Plot');

% Define z-score bins
z_edges = -4:0.1:4;
z_centers = z_edges(1:end-1) + diff(z_edges)/2;

% Compute histogram counts
counts_stim = histcounts(z_STIM, z_edges);
counts_sham = histcounts(z_SHAM, z_edges);
% Convert counts to average per file
counts_stim = counts_stim / num_stim_files;
counts_sham = counts_sham / num_sham_files;

hold on; box on; grid on;
% Plot mirrored horizontal bar plot
barh(z_centers, -counts_sham, 'FaceColor',[0.8 0.3 0.3], 'EdgeColor','none');
barh(z_centers,  counts_stim,  'FaceColor',[0.2 0.6 0.8], 'EdgeColor','none');

xlabel('Number of data points');
ylabel('Z-score (\sigma)');
ylim([-4 4]);

% Add horizontal scale bar (10 peaks)
scale_len = 10;
x_start = -max(counts_sham)*1.1;
y_level = -4.2;
line([x_start, x_start+scale_len], [y_level, y_level], 'Color','k','LineWidth',2);
text(x_start+scale_len/2, y_level-0.1, sprintf('%d peaks', scale_len), ...
     'HorizontalAlignment','center', 'VerticalAlignment','top');

% Add central zero line
line([0 0], ylim, 'Color','k', 'LineWidth',1.5);

legend({'SHAM','STIM'}, 'Location','best');
title('Z-score Distribution Bar Plot');

% Save figures
saveas(fig_zscore, fullfile(outputDir, 'Zscore_BarPlot_Distribution.png'));
saveas(fig_zscore, fullfile(outputDir, 'Zscore_BarPlot_Distribution.fig'));
fprintf('Z-score bar plot saved to: %s\n', fullfile(outputDir, 'Zscore_BarPlot_Distribution.png'));

%% Export summary statistics to console
fprintf('\n=== Analysis Summary ===\n');
if ~isempty(z_STIM) && ~isempty(z_SHAM)
    fprintf('STIM: %d files processed, %d total peaks detected\n', num_stim_files, length(z_STIM));
    fprintf('      Z-score: mean=%.2f, std=%.2f\n', mean(z_STIM), std(z_STIM));
    fprintf('SHAM: %d files processed, %d total peaks detected\n', num_sham_files, length(z_SHAM));
    fprintf('      Z-score: mean=%.2f, std=%.2f\n', mean(z_SHAM), std(z_SHAM));
end

%% Results Summary
fprintf('\nAnalysis complete. Results saved to: %s\n', outputDir);

%% Helper function: Process EEG file and detect peaks
function [peak_times, peak_zscores] = detectPeaksFromFile(filename, params)
    % Load data
    data_struct = load(filename);
    
    % Expect standardized format with 'raw_data' field
    if ~isfield(data_struct, 'raw_data')
        error('Expected field "raw_data" not found in file: %s', filename);
    end
    
    raw_data_matrix = data_struct.raw_data;
    
    % Extract EEG data (first row)
    eeg = double(raw_data_matrix(1,:));
    Fs = 1000; % 1 kHz sampling rate
    
    % Bandpass filter EEG (1-100 Hz)
    [b_eeg, a_eeg] = butter(4, [1 100]/(Fs/2), 'bandpass');
    eeg_filtered = filtfilt(b_eeg, a_eeg, eeg);
    
    % Calculate Moving RMS for detection
    win_samples = round(params.win_sec * Fs);
    if mod(win_samples, 2) == 0, win_samples = win_samples + 1; end
    
    % Calculate squared signal for movmedian
    eeg_sq = eeg_filtered.^2;
    rms_values = sqrt(movmedian(eeg_sq, win_samples, 'omitnan'));
    
    % Set adaptive threshold
    baseline_rms = median(rms_values, 'omitnan');
    threshold = params.k_thresh * baseline_rms;
    
    % Parameters from params
    min_dur_samples = round(params.min_dur_ms / 1000 * Fs);
    refractory_samples = round(params.refractory_ms / 1000 * Fs);
    
    % Find supra-threshold points
    supra_thresh_idx = find(rms_values > threshold);
    
    % Initialize outputs
    peak_times = [];
    peak_intensities = [];
    
    if isempty(supra_thresh_idx)
        peak_times = [];
        peak_zscores = [];
        return; % No peaks found
    end
    
    % Form contiguous runs
    diff_supra = [Inf, diff(supra_thresh_idx)];
    run_starts = find(diff_supra > 1);
    all_run_starts = [1; run_starts(:)];
    
    % Process each run to extract peak times and intensities
    for i = 1:numel(all_run_starts)
        start_idx = all_run_starts(i);
        if i < numel(all_run_starts)
            end_idx = all_run_starts(i+1) - 1;
        else
            end_idx = numel(supra_thresh_idx);
        end
        
        current_run = supra_thresh_idx(start_idx:end_idx);
        
        if numel(current_run) >= min_dur_samples
            % Find peak within this run
            [peak_val, peak_offset] = max(rms_values(current_run));
            peak_idx = current_run(peak_offset);
            
            % Convert to time and store intensity
            peak_time = peak_idx / Fs; % Convert to seconds
            
            % Store results
            peak_times = [peak_times, peak_time];
            peak_intensities = [peak_intensities, peak_val];
        end
    end
    
    % Apply refractory period merging
    if ~isempty(peak_times)
        % Sort by time
        [peak_times, idx] = sort(peak_times);
        peak_intensities = peak_intensities(idx);
        
        % Find peaks that are too close
        i = length(peak_times);
        while i > 1
            if (peak_times(i) - peak_times(i-1)) * Fs < refractory_samples
                % Keep the peak with higher intensity
                if peak_intensities(i) > peak_intensities(i-1)
                    peak_times(i-1) = [];
                    peak_intensities(i-1) = [];
                else
                    peak_times(i) = [];
                    peak_intensities(i) = [];
                end
                i = i - 2;
            else
                i = i - 1;
            end
            
            if i < 1
                break;
            end
        end
    end
    
    % Calculate z-scores based on the mean and std of peaks within this file
    if length(peak_intensities) > 1
        mean_peaks = mean(peak_intensities);
        std_peaks = std(peak_intensities);
        
        if std_peaks > 0
            % Calculate z-scores relative to this file's peak distribution
            peak_zscores = (peak_intensities - mean_peaks) / std_peaks;
        else
            % If std is 0, set z-scores to 0
            peak_zscores = zeros(size(peak_intensities));
        end
    elseif length(peak_intensities) == 1
        % Single peak case - z-score is 0
        peak_zscores = 0;
    else
        % No peaks
        peak_zscores = [];
    end
end