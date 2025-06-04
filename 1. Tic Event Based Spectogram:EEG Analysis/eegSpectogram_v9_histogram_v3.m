%% EEG/EMG Tic Analysis & Visualization (v9 with SNR Peak Analysis)
% Based on user's detailed technical brief for improved peak detection using LF power.

clear; close all; clc;

%% 0) Analysis Parameters (User Tunable)
params.win_sec          = 0.2;   % RMS window (s) for EEG event detection
params.k_thresh         = 3.5;   % Multiplier for EEG RMS threshold (baseline = median(RMS))
params.min_dur_ms       = 50;    % Minimum event duration for EEG (ms)
params.refractory_ms    = 50;   % Refractory period to merge EEG events (ms)

params.LF_band          = [0.5 4];   % Hz, for Low-Frequency power validation
params.LF_win_sec       = 0.25;   % Window (s) for LF power envelope calculation
params.k_LF_thresh      = 3;      % Multiplier for LF power threshold (baseline = median(LF_power))

params.tp_window_ms     = 50;     % Coincidence window for EEG peak and LF power burst (ms)

%% Batch Processing Setup
batch_mode = questdlg('Run in batch mode?', 'Batch Mode Selection', 'Yes', 'No', 'No');
if strcmp(batch_mode, 'Yes')
    % Select multiple files
    [filenames, pathname] = uigetfile('*.mat', 'Select EEG/EMG .mat files', 'MultiSelect', 'on');
    if isequal(filenames, 0)
        disp('User cancelled'); 
        return;
    end
    
    % Convert to cell array if single file selected
    if ~iscell(filenames)
        filenames = {filenames};
    end
    
    % Create output folder for saved figures
    output_folder = fullfile(pathname, 'analyzed_figures');
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
        fprintf('Created output folder: %s\n', output_folder);
    end
    
    fprintf('Batch processing %d files...\n', length(filenames));
else
    % Single file mode
    [filename, pathname] = uigetfile('*.mat','Select EEG/EMG .mat file');
    if isequal(filename,0)
        disp('User cancelled'); 
        return;
    end
    filenames = {filename};
end

%% Process each file
for file_idx = 1:length(filenames)
    % Clear variables from previous iterations
    clear eeg emg eeg_plot_filtered emg_plot_filtered rms_values_eeg EEG_events
    
    current_filename = filenames{file_idx};
    fprintf('\nProcessing file %d/%d: %s\n', file_idx, length(filenames), current_filename);
    
    % Load data
    data_struct = load(fullfile(pathname, current_filename));
    
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
            warning('Could not find suitable EEG data in file: %s. Skipping to next file.', current_filename);
            continue;
        end
    end

    if isempty(raw_data_matrix)
        warning('Loaded data matrix is empty for file: %s. Skipping to next file.', current_filename);
        continue;
    end

    % Initialize variables for proper scope
    has_emg = false;
    emg = [];
    
    eeg = double(raw_data_matrix(1,:));
    if size(raw_data_matrix,1) >= 2
        emg = double(raw_data_matrix(2,:));
        % Additional validation for EMG data
        if ~isempty(emg) && length(emg) == length(eeg) && any(emg ~= 0)
            has_emg = true;
            fprintf('EMG data loaded from row 2 for plotting.\n');
        else
            has_emg = false;
            emg = [];
            fprintf('EMG data invalid or empty. EMG plots will be empty.\n');
        end
    else
        fprintf('No EMG data found in row 2. EMG plots will be empty.\n');
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

    %% Calculate SNR for each event
    fprintf('\n--- Calculating SNR for Each Event ---\n');

    if ~isempty(EEG_events)
        % Add SNR field to EEG_events structure
        [EEG_events.SNR] = deal(0);
        
        % Baseline window length for noise calculation (1 second)
        baseline_window_sec = 1;
        baseline_samples = round(baseline_window_sec * Fs);
        
        for i = 1:numel(EEG_events)
            % Convert times to sample indices
            peak_sample = round(EEG_events(i).peak_time * Fs);
            start_sample = round(EEG_events(i).start_time * Fs);
            end_sample = round(EEG_events(i).end_time * Fs);
            
            % Ensure sample indices are within valid range
            peak_sample = min(max(1, peak_sample), length(eeg_plot_filtered));
            start_sample = min(max(1, start_sample), length(eeg_plot_filtered));
            end_sample = min(max(1, end_sample), length(eeg_plot_filtered));
            
            % Define baseline window for noise calculation
            % Use segment before event start, or if not enough signal, use fixed window elsewhere
            if start_sample > baseline_samples
                % Enough samples before event
                noise_window_indices = (start_sample - baseline_samples):(start_sample - 1);
            else
                % Not enough samples before event, use window after event
                if end_sample + baseline_samples <= length(eeg_plot_filtered)
                    noise_window_indices = (end_sample + 1):(end_sample + baseline_samples);
                else
                    % Use whatever is available
                    available_before = start_sample - 1;
                    available_after = max(0, length(eeg_plot_filtered) - end_sample);
                    if available_before >= available_after && available_before > 0
                        noise_window_indices = 1:available_before;
                    elseif available_after > 0
                        noise_window_indices = (end_sample + 1):length(eeg_plot_filtered);
                    else
                        % Extreme case: use a window from the middle of the recording if event spans whole recording
                        middle_idx = round(length(eeg_plot_filtered)/2);
                        window_size = min(100, length(eeg_plot_filtered));
                        noise_window_indices = max(1, middle_idx - round(window_size/2)):min(middle_idx + round(window_size/2), length(eeg_plot_filtered));
                    end
                end
            end
            
            % Double check indices are valid before accessing array
            noise_window_indices = noise_window_indices(noise_window_indices > 0 & noise_window_indices <= length(eeg_plot_filtered));
            
            % If we ended up with an empty window, use at least something
            if isempty(noise_window_indices)
                noise_window_indices = 1:min(100, length(eeg_plot_filtered));
            end
            
            % Extract the event and noise signals from filtered EEG
            noise_signal = eeg_plot_filtered(noise_window_indices);
            
            % Get peak amplitude at event peak
            peak_amplitude = abs(eeg_plot_filtered(peak_sample));
            
            % Calculate noise standard deviation
            noise_std = std(noise_signal);
            
            % Handle potential zero noise std
            if noise_std == 0
                noise_std = eps; % Use a small value to avoid division by zero
            end
            
            % Calculate SNR in dB
            EEG_events(i).SNR = 20 * log10(peak_amplitude / noise_std);
        end
        
        % Separate SNR values for TP and FP events
        tp_indices = find([EEG_events.isTP]);
        fp_indices = find(~[EEG_events.isTP]);
        
        tp_snr = [EEG_events(tp_indices).SNR];
        fp_snr = [EEG_events(fp_indices).SNR];
        all_snr = [EEG_events.SNR];
        
        % Calculate statistics
        if ~isempty(tp_snr)
            mean_tp_snr = mean(tp_snr);
            fprintf('Mean SNR for True Positive events: %.2f dB\n', mean_tp_snr);
        else
            mean_tp_snr = NaN;
            fprintf('No True Positive events to calculate mean SNR.\n');
        end
        
        if ~isempty(fp_snr)
            mean_fp_snr = mean(fp_snr);
            fprintf('Mean SNR for False Positive events: %.2f dB\n', mean_fp_snr);
        else
            mean_fp_snr = NaN;
            fprintf('No False Positive events to calculate mean SNR.\n');
        end
        
        if ~isempty(all_snr)
            fprintf('Overall mean SNR: %.2f dB\n', mean(all_snr));
        end
    else
        tp_snr = [];
        fp_snr = [];
        all_snr = [];
        fprintf('No events detected to calculate SNR.\n');
    end

    %% Spectrogram Parameters (EEG and EMG)
    winSecSpec   = 0.256; % Window length (s)
    nfftSpec     = round(winSecSpec * Fs);
    noverlapSpec = round(0.9 * nfftSpec); % 90% overlap as per user change
    fMaxEEG_spec = 50;  % Max freq EEG spec (Hz)
    fMaxEMG_spec = 150; % Max freq EMG spec (Hz) - Restored for EMG spectrogram
    caxis_vals_log_uv_sqrt_hz = log10([1e-1 1e3]); % Color axis

    %% Create Combined Figure (2x3 grid with SNR Histogram)
    % Close any existing figures to avoid display issues
    close all;
    
    % Create new figure with proper sizing
    fig_handle = figure('Position', [50, 50, 1600, 800], 'Name', current_filename, ...
                       'NumberTitle', 'off', 'Color', 'white');
    
    % Clear figure and set title
    clf(fig_handle);
    sgtitle_text = sprintf('EEG Event Analysis - File: %s', strrep(current_filename, '_', '\_'));
    sgtitle(sgtitle_text, 'FontWeight','bold', 'FontSize', 14);

    % Create a 2x3 grid layout
    % --- Subplot 1: Filtered EEG + Tics ---
    ax1 = subplot(2,3,1);
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
    ax2 = subplot(2,3,2);
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

    % --- Subplot 3 (NEW): SNR Histogram ---
    ax3 = subplot(2,3,3);
    if ~isempty(EEG_events) && ~isempty(all_snr)
        has_legend_items = false;
        
        hold on;
        if ~isempty(all_snr)
            histogram(all_snr, 'BinWidth', 2, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.3 0.3 0.3], 'DisplayName', 'All Events');
            has_legend_items = true;
        end
        
        if ~isempty(tp_snr)
            histogram(tp_snr, 'BinWidth', 2, 'FaceColor', [0 0.4 0.8], 'EdgeColor', [0 0.2 0.6], 'DisplayName', 'True Positives');
            has_legend_items = true;
        end
        
        if ~isempty(fp_snr)
            histogram(fp_snr, 'BinWidth', 2, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', [0.6 0.1 0.1], 'DisplayName', 'False Positives');
            has_legend_items = true;
        end
        hold off;
        
        xlabel('SNR (dB)', 'FontWeight', 'bold');
        ylabel('Number of Events', 'FontWeight', 'bold');
        title('Distribution of Event Peak SNRs');
        grid on;
        
        if has_legend_items
            legend('show', 'Location', 'best');
        end
        
        % Add statistical annotations
        if ~isempty(tp_snr) || ~isempty(fp_snr)
            yl = ylim;
            xl = xlim;
            text_x = xl(1) + 0.05*(xl(2)-xl(1));
            text_y = yl(2) - 0.1*(yl(2)-yl(1));
            
            if ~isempty(tp_snr) && ~isempty(fp_snr)
                text_stats = sprintf('TP: mean = %.1f dB (n=%d)\nFP: mean = %.1f dB (n=%d)', ...
                    mean_tp_snr, length(tp_snr), mean_fp_snr, length(fp_snr));
            elseif ~isempty(tp_snr)
                text_stats = sprintf('TP: mean = %.1f dB (n=%d)', mean_tp_snr, length(tp_snr));
            elseif ~isempty(fp_snr)
                text_stats = sprintf('FP: mean = %.1f dB (n=%d)', mean_fp_snr, length(fp_snr));
            end
            
            text(text_x, text_y, text_stats, 'FontSize', 9, 'FontWeight', 'bold');
        end
    else
        title('Distribution of Event Peak SNRs (No Events)');
        set(gca, 'Color', [0.95 0.95 0.95]);
        text(0.5, 0.5, 'No events detected', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    % --- Subplot 4: Unfiltered EEG Spectrogram ---
    ax4 = subplot(2,3,4);
    try
        % Validate input signal before spectrogram computation
        if length(eeg) < nfftSpec
            error('EEG signal too short for spectrogram window size');
        end
        
        [fEEG_spec, tEEG_spec, logSpecEEG] = computeLogSpectrogram(eeg, Fs, nfftSpec, noverlapSpec, fMaxEEG_spec);
        
        % Validate spectrogram output
        if isempty(logSpecEEG) || ~all(isfinite(logSpecEEG(:)))
            warning('Invalid spectrogram data, using fallback display');
            logSpecEEG(~isfinite(logSpecEEG)) = min(caxis_vals_log_uv_sqrt_hz);
        end
        
        % Create spectrogram plot
        h_img = imagesc(tEEG_spec, fEEG_spec, logSpecEEG); 
        axis xy;
        colormap(ax4, jet);
        caxis(ax4, caxis_vals_log_uv_sqrt_hz);
        
        % Set proper axis properties
        title('EEG Spectrogram (Unfiltered)', 'FontWeight', 'bold');
        ylabel('Frequency (Hz)', 'FontWeight', 'bold');
        xlabel('Time (s)', 'FontWeight', 'bold');
        
        % Set frequency ticks
        ytick_vals = floor(linspace(0, fMaxEEG_spec, 6));
        set(ax4, 'YTick', ytick_vals, 'TickDir', 'out');
        grid(ax4, 'on'); grid(ax4, 'minor');
        
        fprintf('EEG spectrogram computed successfully (%.1f-%.1f Hz, %.1f-%.1f s)\n', ...
                min(fEEG_spec), max(fEEG_spec), min(tEEG_spec), max(tEEG_spec));
        
    catch ME_eeg_spec
        % Error handling with informative display
        title('EEG Spectrogram - Error', 'Color', 'red', 'FontWeight', 'bold');
        text(0.5, 0.5, sprintf('Error: %s', ME_eeg_spec.message), ...
             'HorizontalAlignment', 'center', 'Color', 'red', 'FontWeight', 'bold');
        set(ax4, 'Color', [0.95 0.95 0.95]);
        axis off;
        warning('Error generating EEG spectrogram: %s', ME_eeg_spec.message);
    end

    % --- Subplot 5: Unfiltered EMG Spectrogram (Restored) ---
    ax5 = subplot(2,3,5);
    if has_emg && ~isempty(emg)
        try
            % Validate input signal before spectrogram computation
            if length(emg) < nfftSpec
                error('EMG signal too short for spectrogram window size');
            end
            
            [fEMG_spec, tEMG_spec, logSpecEMG] = computeLogSpectrogram(emg, Fs, nfftSpec, noverlapSpec, fMaxEMG_spec);
            
            % Validate spectrogram output
            if isempty(logSpecEMG) || ~all(isfinite(logSpecEMG(:)))
                warning('Invalid EMG spectrogram data, using fallback display');
                logSpecEMG(~isfinite(logSpecEMG)) = min(caxis_vals_log_uv_sqrt_hz);
            end
            
            % Create spectrogram plot
            h_img_emg = imagesc(tEMG_spec, fEMG_spec, logSpecEMG); 
            axis xy;
            colormap(ax5, jet);
            caxis(ax5, caxis_vals_log_uv_sqrt_hz);
            
            % Add colorbar for EMG spectrogram
            cbEMG = colorbar(ax5, 'eastoutside');
            ylabel(cbEMG, 'log_{10}(\muV / \surdHz)', 'FontWeight', 'bold');
            
            % Set proper axis properties
            title('EMG Spectrogram (Unfiltered)', 'FontWeight', 'bold');
            ylabel('Frequency (Hz)', 'FontWeight', 'bold');
            xlabel('Time (s)', 'FontWeight', 'bold');
            
            % Set frequency ticks
            ytick_vals = floor(linspace(0, fMaxEMG_spec, 7));
            set(ax5, 'YTick', ytick_vals, 'TickDir', 'out');
            grid(ax5, 'on'); grid(ax5, 'minor');
            
            fprintf('EMG spectrogram computed successfully (%.1f-%.1f Hz, %.1f-%.1f s)\n', ...
                    min(fEMG_spec), max(fEMG_spec), min(tEMG_spec), max(tEMG_spec));
            
        catch ME_emg_spec
            % Error handling with informative display
            title('EMG Spectrogram - Error', 'Color', 'red', 'FontWeight', 'bold');
            text(0.5, 0.5, sprintf('Error: %s', ME_emg_spec.message), ...
                 'HorizontalAlignment', 'center', 'Color', 'red', 'FontWeight', 'bold');
            set(ax5, 'Color', [0.95 0.95 0.95]);
            xlabel('Time (s)', 'FontWeight', 'bold');
            ylabel('Frequency (Hz)', 'FontWeight', 'bold');
            warning('Error generating EMG spectrogram: %s', ME_emg_spec.message);
        end
    else
        % No EMG data available
        title('EMG Spectrogram (No Data)', 'FontWeight', 'bold');
        text(0.5, 0.5, 'No EMG Data Available', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'FontSize', 12);
        set(ax5, 'Color', [0.95 0.95 0.95]);
        xlabel('Time (s)', 'FontWeight', 'bold');
        ylabel('Frequency (Hz)', 'FontWeight', 'bold');
        axis([0 1 0 1]); % Set standard axis limits for empty plot
    end
    
    % --- Subplot 6: EEG Band Power Distribution ---
    subplot(2,3,6);
    ax6 = gca; % Get current axis handle

    % Initialize variables to avoid undefined references in catch block
    band_names = {};
    band_power = [];
    EEG_band_ratios = struct();
    
    try
        % Calculate EEG frequency bands similar to FFT_psd_histogram_v3.m
        % Define frequency bands
        bands = {'Delta', [0.5,4];
                 'Theta', [4,7]; 
                 'Alpha', [8,12];
                 'Beta',  [12,30];
                 'Gamma', [30,100]};
        
        num_bands = size(bands, 1);
        band_names = bands(:,1);
        
        % Calculate power for each band using bandpower function
        band_power = zeros(1, num_bands);
        for b = 1:num_bands
            freq_range = bands{b,2};
            try
                band_power(b) = bandpower(eeg, Fs, freq_range);
            catch band_err
                % Fallback if bandpower function is not available
                fprintf('Warning: bandpower function not available, using manual FFT calculation\n');
                % Manual calculation using FFT
                L = length(eeg);
                Y = fft(eeg);
                P2 = abs(Y/L);
                P1 = P2(1:floor(L/2)+1);
                P1(2:end-1) = 2*P1(2:end-1);
                freq = Fs*(0:(floor(L/2)))/L;
                
                % Find indices for this band
                band_idx = (freq >= freq_range(1)) & (freq <= freq_range(2));
                band_power(b) = sum(P1(band_idx).^2);
            end
        end
        
        % Calculate percentages
        total_power = sum(band_power);
        if total_power > 0
            band_percentages = (band_power / total_power) * 100;
        else
            band_percentages = zeros(1, num_bands);
        end
        
        % Create side-by-side bar chart with custom colors
        band_colors = [0.2 0.6 0.8;   % Delta - blue
                      0.8 0.2 0.2;   % Theta - red
                      0.2 0.8 0.2;   % Alpha - green
                      0.8 0.8 0.2;   % Beta - yellow
                      0.8 0.2 0.8];  % Gamma - purple
        
        % Create bar chart showing percentages
        h = bar(band_percentages, 'FaceColor', 'flat');
        
        % Apply colors
        for i = 1:length(bands)
            h.CData(i,:) = band_colors(i,:);
        end
        
        % Set axis properties
        set(ax6, 'XTick', 1:num_bands, 'XTickLabel', band_names);
        xtickangle(45);
        ylim([0 max(band_percentages)*1.2]); % Give some headroom
        ylabel('Power (%)');
        title('EEG Band Power Distribution');
        grid on;
        
        % Add percentages as text on top of each bar
        for i = 1:num_bands
            % Use text function with explicit parent to avoid deletion issues
            text(i, band_percentages(i) + 1, sprintf('%.1f%%', band_percentages(i)), ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        
        % Add ratio calculations (common in EEG analysis)
        delta_idx = find(strcmp(band_names, 'Delta'));
        theta_idx = find(strcmp(band_names, 'Theta'));
        alpha_idx = find(strcmp(band_names, 'Alpha'));
        beta_idx = find(strcmp(band_names, 'Beta'));
        gamma_idx = find(strcmp(band_names, 'Gamma'));
        
        % Calculate EEG ratios (stored in variables but not displayed)
        if band_power(beta_idx) > 0
            theta_beta_ratio = band_power(theta_idx) / band_power(beta_idx);
        else
            theta_beta_ratio = NaN;
        end
        
        if band_power(delta_idx) > 0
            theta_delta_ratio = band_power(theta_idx) / band_power(delta_idx);
        else
            theta_delta_ratio = NaN;
        end
        
        % Store ratios in EEG_band_ratios struct for potential future use
        EEG_band_ratios.theta_beta = theta_beta_ratio;
        EEG_band_ratios.theta_delta = theta_delta_ratio;
        
    catch ME
        % Handle any errors in this subplot creation
        warning('Error creating EEG band power distribution plot: %s', ME.message);
        title('EEG Band Power Distribution - Error');
        text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
             'HorizontalAlignment', 'center', 'Color', 'red');
        axis off;
    end
    
    % --- Link Axes and Final Adjustments ---
    % Collect axes handles that should be time-linked (exclude histogram and band power plots)
    axes_to_link = [];
    
    % Only include valid axes that have time-based data
    if exist('ax1', 'var') && isvalid(ax1) && isgraphics(ax1)
        axes_to_link = [axes_to_link, ax1];
    end
    if exist('ax2', 'var') && isvalid(ax2) && isgraphics(ax2)
        axes_to_link = [axes_to_link, ax2];
    end
    if exist('ax4', 'var') && isvalid(ax4) && isgraphics(ax4)
        axes_to_link = [axes_to_link, ax4];
    end
    if exist('ax5', 'var') && isvalid(ax5) && isgraphics(ax5) && has_emg
        axes_to_link = [axes_to_link, ax5];
    end
    
    % Link axes only if we have valid axes to link
    if length(axes_to_link) > 1
        try
            linkaxes(axes_to_link, 'x');
            full_time_range = [t(1), t(end)]; 
            set(axes_to_link, 'XLim', full_time_range);
            fprintf('Successfully linked %d time-based axes\n', length(axes_to_link));
        catch link_err
            warning('Could not link axes: %s', link_err.message);
        end
    end

    % Ensure figure is properly sized and positioned
    try
        set(fig_handle, 'Position', [50, 50, 1600, 800]);
        
        % Improve subplot spacing
        sgtitle(sgtitle_text, 'FontWeight','bold', 'FontSize', 14);
        
        % Force figure update and bring to front
        drawnow; 
        figure(fig_handle); % Bring figure to front
        
        fprintf('Figure display completed successfully\n');
    catch fig_err
        warning('Could not finalize figure display: %s', fig_err.message);
    end
    
    % Output event counts and rate to console
    total_duration_min = t(end) / 60;
    if total_duration_min > 0
        event_rate_per_min = num_artifacts_eeg / total_duration_min;
    else
        event_rate_per_min = 0;
    end
    
    fprintf('Final EEG Events Summary:\nTotal Detected: %d\nTrue Positives (EEG && LF): %d\nFalse Positives (EEG only): %d\nEvent Rate: %.2f events/min\n', ...
            num_artifacts_eeg, tp_count, fp_count, event_rate_per_min);

    % Save figure if in batch mode
    if strcmp(batch_mode, 'Yes')
        % Create filename for the saved figure
        [~, filename_noext, ~] = fileparts(current_filename);
        save_filename = fullfile(output_folder, [filename_noext '_analyzed.fig']);
        savefig(gcf, save_filename);
        fprintf('Saved figure to: %s\n', save_filename);
    end
    
    % Pause to allow user to view the figure if there are multiple files
    if strcmp(batch_mode, 'Yes') && file_idx < length(filenames)
        pause(0.5); % Brief pause between files
    end
end

% Note: computeLogSpectrogram has been moved to a separate file
