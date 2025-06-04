%% EEG/EMG Tic Analysis & Visualization (v8 re-engineered with LF Validation)
% Based on user's detailed technical brief for improved peak detection using LF power.

clear; close all; clc;

%% 0) Analysis Parameters (User Tunable)
params.win_sec          = 0.5;   % RMS window (s) for EEG event detection
params.k_thresh         = 3.5;      % Multiplier for EEG RMS threshold (baseline = median(RMS))
params.min_dur_ms       = 50;     % Minimum event duration for EEG (ms)
params.refractory_ms    = 100;    % Refractory period to merge EEG events (ms)

params.LF_band          = [0.5 4];   % Hz, for Low-Frequency power validation
params.LF_win_sec       = 0.25;   % Window (s) for LF power envelope calculation
params.k_LF_thresh      = 3;      % Multiplier for LF power threshold (baseline = median(LF_power))

params.tp_window_ms     = 50;     % Coincidence window for EEG peak and LF power burst (ms)

%% Load EEG Data (EMG loading restored for plotting if available)
[filename, pathname] = uigetfile('*.mat','Select EEG/EMG .mat file');
if isequal(filename,0), disp('User cancelled'); return; end
data_struct = load(fullfile(pathname, filename));

% Flexible data loading (assuming EEG=row1)
raw_data_matrix = [];
if isfield(data_struct, 'raw_data')
    raw_data_matrix = data_struct.raw_data;
    fprintf('Data loaded from field: raw_data\n');
elseif isfield(data_struct, 'data')
    raw_data_matrix = data_struct.data;
    fprintf('Data loaded from field: data\n');
else
    fields = fieldnames(data_struct);
    found_data = false;
    for i = 1:length(fields)
        current_field_data = data_struct.(fields{i});
        if isnumeric(current_field_data) && ismatrix(current_field_data) && size(current_field_data,1) >= 1
            raw_data_matrix = current_field_data;
            fprintf('Data loaded from field: %s\n', fields{i});
            found_data = true;
            break;
        end
    end
    if ~found_data
        errordlg('Could not find suitable EEG data in the .mat file.', 'Data Loading Error');
        return;
    end
end

if isempty(raw_data_matrix)
    errordlg('Loaded data matrix is empty.', 'Data Error');
    return;
end

eeg = double(raw_data_matrix(1,:));
if size(raw_data_matrix,1) >= 2
    emg = double(raw_data_matrix(2,:));
    fprintf('EMG data loaded from row 2 for plotting.\n');
    has_emg = true;
else
    emg = []; 
    fprintf('No EMG data found in row 2. EMG plots will be empty.\n');
    has_emg = false;
end

Fs = 1000; % Assuming 1 kHz sampling rate
t = (0:length(eeg)-1)/Fs;

%% Bandpass filter EEG (1-100 Hz for plotting and initial detection base)
[b_eeg_plot, a_eeg_plot] = butter(4, [1 100]/(Fs/2), 'bandpass');
eeg_plot_filtered = filtfilt(b_eeg_plot, a_eeg_plot, eeg);
eeg_detection_filtered = eeg_plot_filtered;

if has_emg
    emg_plot_filtered = filtfilt(b_eeg_plot, a_eeg_plot, emg); % Filter EMG for plotting with 1-100Hz
    fprintf('Applied 1-100Hz filter to EEG and EMG (for plotting).\n');
else
    fprintf('Applied 1-100Hz filter to EEG (for plotting and detection base).\n');
end

%% Calculate Moving RMS for EEG Event Detection (Median-based)
% Uses eeg_detection_filtered (1-100Hz BPF EEG)
win_samples_eeg = round(params.win_sec * Fs);
if mod(win_samples_eeg, 2) == 0, win_samples_eeg = win_samples_eeg + 1; end % Ensure odd for movmedian

% Calculate squared signal for movmedian
eeg_sq = eeg_detection_filtered.^2;
rms_values_eeg = sqrt(movmedian(eeg_sq, win_samples_eeg, 'omitnan')); % Use omitnan for robustness at edges

%% Set Adaptive RMS Threshold for EEG
baseline_eeg_rms = median(rms_values_eeg, 'omitnan');
threshold_eeg = params.k_thresh * baseline_eeg_rms;

%% (NEW) EEG Event Extraction (Two-Pass: Duration Filter then Refractory Merge)

% --- Parameters from params struct ---
min_dur_samples_eeg = round(params.min_dur_ms / 1000 * Fs);
refractory_samples_eeg = round(params.refractory_ms / 1000 * Fs);

EEG_events = struct('start_time', {}, 'peak_time', {}, 'end_time', {}, ...
                    'duration_ms', {}, 'peak_rms_value', {}, 'isTP', {});

% 1. Find supra-threshold points
supra_thresh_idx_eeg = find(rms_values_eeg > threshold_eeg);

if ~isempty(supra_thresh_idx_eeg)
    % 2. Form contiguous runs (gap <= 1 sample for initial run formation)
    diff_supra = [Inf, diff(supra_thresh_idx_eeg)];
    run_starts_in_supra_idx = find(diff_supra > 1);
    all_run_starts_eeg = [1; run_starts_in_supra_idx(:)]; % Ensure column and starts with 1st supra_thresh_idx
    
    candidate_events = []; % Store events passing duration filter
    
    for i = 1:numel(all_run_starts_eeg)
        start_of_run_idx_in_supra = all_run_starts_eeg(i);
        if i < numel(all_run_starts_eeg)
            end_of_run_idx_in_supra = all_run_starts_eeg(i+1) - 1;
        else
            end_of_run_idx_in_supra = numel(supra_thresh_idx_eeg);
        end
        
        current_run_indices = supra_thresh_idx_eeg(start_of_run_idx_in_supra : end_of_run_idx_in_supra);
        
        if numel(current_run_indices) >= min_dur_samples_eeg
            % Event is long enough, find peak RMS within this run
            [peak_val, peak_offset_in_run] = max(rms_values_eeg(current_run_indices));
            peak_idx_eeg = current_run_indices(peak_offset_in_run);
            
            candidate_event.start_idx = current_run_indices(1);
            candidate_event.peak_idx  = peak_idx_eeg;
            candidate_event.end_idx   = current_run_indices(end);
            candidate_event.peak_rms_value = peak_val;
            candidate_events = [candidate_events, candidate_event];
        end
    end
    
    % 3. Merge candidate events based on refractory period (inter-PEAK interval)
    if ~isempty(candidate_events)
        % Sort candidates by peak time to make merging easier
        [~, sort_order] = sort([candidate_events.peak_idx]);
        candidate_events = candidate_events(sort_order);
        
        merged_event_idx = 1;
        EEG_events_final(merged_event_idx) = candidate_events(1); % Start with the first event
        
        for i = 2:numel(candidate_events)
            current_event = candidate_events(i);
            last_merged_event = EEG_events_final(merged_event_idx);
            
            inter_peak_interval_samples = current_event.peak_idx - last_merged_event.peak_idx;
            
            if inter_peak_interval_samples < refractory_samples_eeg
                % Merge: extend the last merged event to include the current one
                EEG_events_final(merged_event_idx).end_idx = max(last_merged_event.end_idx, current_event.end_idx);
                % Update peak if current event has a higher peak RMS
                if current_event.peak_rms_value > last_merged_event.peak_rms_value
                    EEG_events_final(merged_event_idx).peak_idx = current_event.peak_idx;
                    EEG_events_final(merged_event_idx).peak_rms_value = current_event.peak_rms_value;
                end
                % Start time remains from the earlier of the two (already sorted by peak, so start of last_merged)
            else
                % No merge, this is a new distinct event
                merged_event_idx = merged_event_idx + 1;
                EEG_events_final(merged_event_idx) = current_event;
            end
        end
        
        % Convert final merged event indices to times and populate EEG_events struct
        for i = 1:numel(EEG_events_final)
            EEG_events(i).start_time   = EEG_events_final(i).start_idx / Fs;
            EEG_events(i).peak_time    = EEG_events_final(i).peak_idx / Fs;
            EEG_events(i).end_time     = EEG_events_final(i).end_idx / Fs;
            EEG_events(i).duration_ms  = (EEG_events_final(i).end_idx - EEG_events_final(i).start_idx + 1) / Fs * 1000;
            EEG_events(i).peak_rms_value = EEG_events_final(i).peak_rms_value;
            EEG_events(i).isTP         = false; % Default to False Positive
        end
    end
end

num_artifacts_eeg = numel(EEG_events);
fprintf('Initial EEG events detected (after duration & refractory): %d\n', num_artifacts_eeg);
fprintf('EEG RMS threshold: %.2f uV, Baseline EEG RMS: %.2f uV\n', threshold_eeg, baseline_eeg_rms);

%% (NEW) Low-Frequency (LF) Power Module for True Positive Validation
fprintf('\n--- Calculating Low-Frequency Power for Validation ---\n');

% 1. Band-pass eeg_detection_filtered to LF_band
[b_LF, a_LF] = butter(2, [0.5 10]/(Fs/2), 'bandpass');

% Apply padding for filtfilt
pad_samples = 500;
if length(eeg_detection_filtered) > 2 * pad_samples % Ensure signal is long enough for padding
    % Corrected padding for row vector eeg_detection_filtered
    start_padding_segment = eeg_detection_filtered(1:pad_samples);
    end_padding_segment = eeg_detection_filtered(end-pad_samples+1:end);
    eeg_padded = [fliplr(start_padding_segment), eeg_detection_filtered, fliplr(end_padding_segment)];
    
    eeg_LF_padded = filtfilt(b_LF, a_LF, eeg_padded);
    eeg_LF = eeg_LF_padded(pad_samples+1:end-pad_samples);   % remove the padding
else
    warning('Signal too short for specified padding of %d samples. Applying filter without padding.', pad_samples);
    eeg_LF = filtfilt(b_LF, a_LF, eeg_detection_filtered);
end

% Insert sanity guard
if all(~isfinite(eeg_LF))
      error('LF band-pass produced only NaNs/Infs – adjust filter.');
end

% Diagnostic print for eeg_LF
disp('eeg_LF diagnostics [min, max, sum(isnan)]:');
disp([min(eeg_LF) max(eeg_LF) sum(isnan(eeg_LF))]);

% 2. Instantaneous Power Envelope
win_samples_LF = round(params.LF_win_sec*Fs);
if mod(win_samples_LF, 2) == 0, win_samples_LF = win_samples_LF+1; end % Ensure odd for movmean/movmedian
LF_power = movmean( eeg_LF.^2, win_samples_LF, 'omitnan'); % As per brief, or movmedian

% 3. Adaptive Threshold for LF Power
finite_LF = LF_power(isfinite(LF_power));
if isempty(finite_LF) % Handle case where all LF_power is NaN after filtering
    warning('All LF_power values are NaN after filtering and movmean. LF_baseline will be NaN.');
    LF_baseline = NaN;
else
    LF_baseline = median(finite_LF);
end
LF_thr   = params.k_LF_thresh * LF_baseline;

% 4. Binary Mask of LF Power Bursts
LF_burst_idx_mask = LF_power > LF_thr;

fprintf('LF Power threshold: %.2f (based on k=%.1f * median LF power=%.2f)\n', LF_thr, params.k_LF_thresh, LF_baseline);

%% (NEW) Label True Tics based on LF Power Coincidence
if ~isempty(EEG_events)
    tp_samples = round(params.tp_window_ms / 1000 * Fs);
    num_eeg_events_for_tp_check = numel(EEG_events);

    for i = 1:num_eeg_events_for_tp_check
        peak_sample_eeg = round(EEG_events(i).peak_time * Fs);
        
        % Define window around EEG peak for checking LF burst coincidence
        start_check_idx = max(1, peak_sample_eeg - tp_samples);
        end_check_idx   = min(length(LF_burst_idx_mask), peak_sample_eeg + tp_samples);
        
        window_indices = start_check_idx:end_check_idx;
        
        if any(LF_burst_idx_mask(window_indices))
            EEG_events(i).isTP = true;
        end
    end
    fprintf('True Positive labeling complete based on LF power coincidence.\n');
else
    fprintf('No EEG events detected to label for True Positives.\n');
end

% Calculate final TP/FP counts
tp_count = sum([EEG_events.isTP]);
fp_count = num_artifacts_eeg - tp_count;

%% Spectrogram Parameters (EEG and EMG)
winSecSpec   = 0.256; % Window length (s)
nfftSpec     = round(winSecSpec * Fs);
noverlapSpec = round(0.9 * nfftSpec); % 90% overlap as per user change
fMaxEEG_spec = 50;  % Max freq EEG spec (Hz)
fMaxEMG_spec = 150; % Max freq EMG spec (Hz) - Restored for EMG spectrogram
caxis_vals_log_uv_sqrt_hz = log10([1e-1 1e3]); % Color axis

%% Create Combined Figure (2x2)
figure(1);
clf;
sgtitle_text = sprintf('EEG Event Analysis (LF Validated) & EMG View - File: %s', strrep(filename, '_', '\_'));
sgtitle(sgtitle_text, 'FontWeight','bold');

% --- Subplot 1: Filtered EEG + Tics ---
ax1 = subplot(2,2,1);
plot(t, eeg_plot_filtered, 'k'); hold on;

% Plot all EEG events (Red dots)
if ~isempty(EEG_events)
    all_peak_times = [EEG_events.peak_time];
    plot(all_peak_times, zeros(size(all_peak_times)), 'ro', ...
        'MarkerFaceColor','r', 'MarkerSize', 6, 'LineWidth',1, 'DisplayName', 'EEG Event (All)');
end

% Plot True Positive EEG events (Blue outline)
EEG_TP_indices = find([EEG_events.isTP]);
if ~isempty(EEG_TP_indices)
    tp_peak_times = [EEG_events(EEG_TP_indices).peak_time];
    plot(tp_peak_times, zeros(size(tp_peak_times)), 'bo', ...
        'MarkerFaceColor','none', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'EEG Event (LF Validated)');
end

hold off;
axis tight; grid on;
title('Filtered EEG (1-100 Hz) with Detected Events');
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
legend('show', 'Location','northwest');

% Annotate event count on this subplot
xl = xlim; yl = ylim;
text_str = sprintf('EEG Events: %d (TP: %d, FP: %d)', num_artifacts_eeg, tp_count, fp_count);
text(xl(1) + 0.02*(xl(2)-xl(1)), yl(2) - 0.1*(yl(2)-yl(1)), text_str, 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

% --- Subplot 2: Filtered EMG (Restored) ---
ax2 = subplot(2,2,2);
if has_emg
    plot(t, emg_plot_filtered, 'Color', [0.4 0.4 0.4]);
    title('Filtered EMG (1-100 Hz)');
    ylabel('Amplitude (\muV)');
else
    plot(t, nan(size(t)),'w.'); % Plot invisible point to keep axes structure
    title('Filtered EMG (No Data)');
    ylabel('Amplitude (\muV)');
    set(gca, 'Color', [0.95 0.95 0.95], 'XTickLabel', {}, 'YTickLabel', {});
end
axis tight; grid on;
xlabel('Time (s)');

% --- Subplot 3: Unfiltered EEG Spectrogram ---
ax3 = subplot(2,2,3);
try
    [fEEG_spec, tEEG_spec, logSpecEEG] = computeLogSpectrogram(eeg, Fs, nfftSpec, noverlapSpec, fMaxEEG_spec); % Use original UNFILTERED eeg
    imagesc(tEEG_spec, fEEG_spec, logSpecEEG); axis xy;
    colormap(ax3, jet);
    caxis(ax3, caxis_vals_log_uv_sqrt_hz);
    % Colorbar for EEG spec removed as per user request
    title('EEG Spectrogram (Unfiltered)');
    ylabel('Frequency (Hz)');
    set(ax3, 'YTick', floor(linspace(0,fMaxEEG_spec,6)), 'TickDir','out');
catch ME_eeg_spec
    title('EEG Spectrogram - Error');
    text(0.5,0.5, sprintf('Error: %s', ME_eeg_spec.message), 'HorizontalAlignment', 'center');
    warning('Error generating EEG spectrogram: %s', ME_eeg_spec.message);
end
xlabel('Time (s)');

% --- Subplot 4: Unfiltered EMG Spectrogram (Restored) ---
ax4 = subplot(2,2,4);
if has_emg
    try
        [fEMG_spec, tEMG_spec, logSpecEMG] = computeLogSpectrogram(emg, Fs, nfftSpec, noverlapSpec, fMaxEMG_spec); % Use original UNFILTERED emg
        imagesc(tEMG_spec, fEMG_spec, logSpecEMG); axis xy;
        colormap(ax4, jet);
        caxis(ax4, caxis_vals_log_uv_sqrt_hz);
        cbEMG = colorbar(ax4); % Restore individual colorbar for EMG spectrogram
        ylabel(cbEMG,'log_{10}(\muV / \surdHz)');
        title('EMG Spectrogram (Unfiltered)');
        ylabel('Frequency (Hz)');
        set(ax4, 'YTick', floor(linspace(0,fMaxEMG_spec,7)), 'TickDir','out');
    catch ME_emg_spec
        title('EMG Spectrogram - Error');
        text(0.5,0.5, sprintf('Error: %s', ME_emg_spec.message), 'HorizontalAlignment', 'center');
        warning('Error generating EMG spectrogram: %s', ME_emg_spec.message);
    end
else
    title('EMG Spectrogram (No Data)');
    set(gca, 'Color', [0.95 0.95 0.95], 'XTickLabel', [], 'YTickLabel', []);
    axis off; 
end
xlabel('Time (s)');

% --- Link Axes and Final Adjustments ---
all_axes = [ax1, ax3]; % Start with essential EEG axes
if has_emg % Only add EMG axes to linking if they contain data
    all_axes = [all_axes, ax2, ax4];
else % If no EMG, ax2 and ax4 might not be fully interactive, link what makes sense
    all_axes = [all_axes, ax2]; % ax2 is created even if blank, ax4 might be fully off
    if isgraphics(ax4) && isvalid(ax4) % if ax4 was attempted but failed or is truly blank
         % all_axes = [all_axes, ax4]; % Decide if linking a blank ax4 is useful
    end
end
% Clean up all_axes to only include valid graphics handles before linking
valid_axes_to_link = all_axes(isgraphics(all_axes) & arrayfun(@isvalid, all_axes));
if ~isempty(valid_axes_to_link)
    linkaxes(valid_axes_to_link, 'x');
    full_time_range = [t(1), t(end)]; 
    set(valid_axes_to_link, 'XLim', full_time_range);
end

set(gcf, 'Position', [50, 50, 1400, 800]); % Restore larger figure size
disp('Combined analysis figure generated.');

% Output event counts and rate to console
total_duration_min = t(end) / 60;
if total_duration_min > 0
    event_rate_per_min = num_artifacts_eeg / total_duration_min;
else
    event_rate_per_min = 0;
end
fprintf('Final EEG Events Summary:\nTotal Detected: %d\nTrue Positives (EEG && LF): %d\nFalse Positives (EEG only): %d\nEvent Rate: %.2f events/min\n', ...
        num_artifacts_eeg, tp_count, fp_count, event_rate_per_min);


%% Helper Function Definition for Spectrogram (from v6)
function [fOut, tOut, logSpec] = computeLogSpectrogram(sig, Fs_in, NFFT_val, NOverlap_val, fCut)
    % Use hamming window for better spectral estimation
    window_vec = hamming(NFFT_val);
    
    % MATLAB's spectrogram with 'psd' gives Power Spectral Density
    [~, fAll_psd, tOut, P_psd] = spectrogram(double(sig), window_vec, NOverlap_val, NFFT_val, Fs_in, 'psd');
    
    % Select frequencies up to fCut
    keep_freqs = fAll_psd <= fCut;
    fOut = fAll_psd(keep_freqs);
    P_subset = P_psd(keep_freqs,:);
    
    % Floor PSD values to prevent log(0) or log(negative) for ASD calculation
    % Smallest ASD ~1e-6 µV/√Hz => smallest PSD ~1e-12 (µV)^2/Hz
    min_psd_val = 1e-12;
    P_subset(P_subset < min_psd_val) = min_psd_val;
    
    % Compute log10(ASD)
    % ASD = sqrt(PSD); units: µV/√Hz if input sig is in µV.
    logSpec = log10(sqrt(P_subset));
end
