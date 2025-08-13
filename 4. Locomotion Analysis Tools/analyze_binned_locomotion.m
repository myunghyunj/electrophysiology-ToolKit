% analyze_binned_locomotion.m
% Analyze locomotion data from head tracking results
% Creates binned distance bar plots for SHAM and STIM groups

%% USER PARAMETERS
bin_size_min = input('Enter bin size in minutes (e.g., 1, 5, 10): ');
mm_per_pix   = input('Enter mm per pixel (e.g., 0.533): ');

if isempty(bin_size_min) || bin_size_min<=0
    error('Invalid bin size. Must be positive.');
end
if isempty(mm_per_pix) || mm_per_pix<=0
    error('Invalid mm‑per‑pixel value.');
end

fprintf('\n=== Binned Locomotion Analysis ===\n');
fprintf('Bin size: %d minutes\n', bin_size_min);

%% GATHER OUTPUT FILES
% Look for head tracking output files
files = dir('*_head_tracking_3fps.xlsx');
if isempty(files)
    % Check in subdirectory
    files = dir('head_tracking_results/*_head_tracking_3fps.xlsx');
    if isempty(files)
        error('No head tracking output files found. Run dlc_head_tracking_resampled.m first.');
    end
end

% Separate SHAM and STIM files
sham_files = {};
stim_files = {};

for i = 1:numel(files)
    if contains(upper(files(i).name), 'SHAM')
        sham_files{end+1} = fullfile(files(i).folder, files(i).name);
    elseif contains(upper(files(i).name), 'STIM')
        stim_files{end+1} = fullfile(files(i).folder, files(i).name);
    end
end

fprintf('Found %d SHAM files and %d STIM files\n\n', length(sham_files), length(stim_files));

%% PROCESS FILES BY GROUP
[sham_binned_data, sham_max_bins] = process_group(sham_files, bin_size_min, mm_per_pix,'SHAM');
[stim_binned_data, stim_max_bins] = process_group(stim_files, bin_size_min, mm_per_pix,'STIM');

% Determine the maximum number of bins across both groups
max_bins = max(sham_max_bins, stim_max_bins);

fprintf('Maximum bins: SHAM=%d, STIM=%d, Overall=%d\n', sham_max_bins, stim_max_bins, max_bins);

%% CREATE VISUALIZATION
fig = figure('Position', [100, 100, 1400, 600], 'Name', 'Binned Locomotion Analysis');

% SHAM subplot (left)
subplot(1, 2, 1);
create_group_plot(sham_binned_data, bin_size_min, max_bins, 'SHAM', [0.3, 0.5, 0.9]);

% STIM subplot (right)
subplot(1, 2, 2);
create_group_plot(stim_binned_data, bin_size_min, max_bins, 'STIM', [0.9, 0.3, 0.3]);

% Save figure
saveas(fig, sprintf('binned_locomotion_%dmin.png', bin_size_min));
fprintf('\nFigure saved as: binned_locomotion_%dmin.png\n', bin_size_min);

%% STATISTICAL SUMMARY
fprintf('\n=== Statistical Summary ===\n');
print_group_stats(sham_binned_data, 'SHAM', bin_size_min);
print_group_stats(stim_binned_data, 'STIM', bin_size_min);

% Export summary data
export_summary(sham_binned_data, stim_binned_data, bin_size_min, max_bins);

%% HELPER FUNCTIONS
function [binned_data, max_bins] = process_group(file_list, bin_size_min, mm_per_pix, group_name)
    fprintf('Processing %s files...\n', group_name);
    binned_data = {};
    max_bins = 0;
    
    for i = 1:length(file_list)
        % Read data from Excel file
        data = readtable(file_list{i});
        
        % Extract filename for display
        [~, basename, ~] = fileparts(file_list{i});
        
        % Convert time to minutes
        time_min = data.Time_sec / 60;
        distances = data.Distance_pixels * mm_per_pix;   % ← now millimetres
        
        % Determine number of bins
        total_time_min = max(time_min);
        n_bins = ceil(total_time_min / bin_size_min);
        
        % Safety check: ensure n_bins is at least 1
        if n_bins < 1
            n_bins = 1;
            warning('File %s has very short duration, setting n_bins to 1', basename);
        end
        
        max_bins = max(max_bins, n_bins);
        
        % Bin the data
        binned_distances = zeros(n_bins, 1);
        
        for b = 1:n_bins
            bin_start = (b - 1) * bin_size_min;
            bin_end = b * bin_size_min;
            
            % Find indices within this bin
            bin_idx = find(time_min >= bin_start & time_min < bin_end);
            
            if ~isempty(bin_idx)
                binned_distances(b) = sum(distances(bin_idx));
            end
        end
        
        % Store results
        binned_data{i}.filename = basename;
        binned_data{i}.distances = binned_distances;
        binned_data{i}.n_bins = n_bins;
        
        % Safety check: ensure distances array is not empty
        if isempty(binned_distances) || all(binned_distances == 0)
            warning('File %s has no valid distance data', basename);
        end
        
        fprintf('  %s: %d bins\n', basename, n_bins);
    end
end

function create_group_plot(binned_data, bin_size_min, max_bins, group_name, color)
    if isempty(binned_data)
        text(0.5, 0.5, sprintf('No %s data', group_name), ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
        return;
    end
    
    % Calculate mean and SEM for each bin
    all_distances = NaN(length(binned_data), max_bins);
    
    for i = 1:length(binned_data)
        n_bins = binned_data{i}.n_bins;
        all_distances(i, 1:n_bins) = binned_data{i}.distances;
    end
    
    % Calculate statistics
    mean_distances = nanmean(all_distances, 1);
    sample_counts  = sum(~isnan(all_distances),1);       % 1×max_bins
    sample_counts(sample_counts==0) = 1;                 % guard
    sem_distances  = nanstd(all_distances,0,1) ./ sqrt(sample_counts);
    
    % Create bar plot
    x_labels = (1:max_bins) * bin_size_min;
    bar_handles = bar(x_labels, mean_distances, 'FaceColor', color, 'EdgeColor', 'k');
    
    % Add error bars
    hold on;
    errorbar(x_labels, mean_distances, sem_distances, 'k.', 'LineWidth', 1.5);
    
    % Customize plot
    xlabel('Time (minutes)');
    ylabel('Distance (mm)');
    title(sprintf('%s Group (n=%d)', group_name, length(binned_data)));
    grid on;
    
    % Set x-axis ticks
    xticks(x_labels);
    xticklabels(arrayfun(@(x) sprintf('%d-%d', x-bin_size_min, x), x_labels, 'UniformOutput', false));
    xtickangle(45);
    
    % Add individual data points
    for i = 1:length(binned_data)
        n_bins = binned_data{i}.n_bins;
        scatter(x_labels(1:n_bins), binned_data{i}.distances, 20, ...
            'k', 'filled', 'MarkerFaceAlpha', 0.3);
    end
    
    % Add sample size annotation
    for b = 1:max_bins
        n_samples = sum(~isnan(all_distances(:, b)));
        if n_samples > 0
            text(x_labels(b), max(mean_distances(b) + sem_distances(b), 0) * 1.05, ...
                sprintf('n=%d', n_samples), 'HorizontalAlignment', 'center', ...
                'FontSize', 8);
        end
    end
end

function print_group_stats(binned_data, group_name, bin_size_min)
    if isempty(binned_data)
        fprintf('\n%s Group: No data\n', group_name);  return
    end

    fprintf('\n%s Group Statistics:\n', group_name);
    fprintf('  Number of animals: %d\n', numel(binned_data));

    % total distance per animal (mm)
    total_distances = cellfun(@(x) sum(x.distances), binned_data);

    fprintf('  Total distance (mean ± SEM): %.1f ± %.1f mm\n', ...
            mean(total_distances), ...
            std(total_distances)/sqrt(numel(total_distances)));
    fprintf('  Total distance range: %.1f – %.1f mm\n', ...
            min(total_distances), max(total_distances));

    % mean distance per bin (mm)
    bin_means = cellfun(@(x) mean(x.distances), binned_data);
    fprintf('  Average distance per %d‑min bin: %.1f ± %.1f mm\n\n', ...
            bin_size_min, mean(bin_means), ...
            std(bin_means)/sqrt(numel(bin_means)));
end


function export_summary(sham_data, stim_data, bin_size_min, max_bins)
    % Create summary table
    summary = table();
    
    % Add bin information
    bin_labels = cell(max_bins, 1);
    for b = 1:max_bins
        bin_labels{b} = sprintf('%d-%d min', (b-1)*bin_size_min, b*bin_size_min);
    end
    summary.Bin = bin_labels;
    
    % Calculate group means
    if ~isempty(sham_data)
        sham_all = NaN(length(sham_data), max_bins);
        for i = 1:length(sham_data)
            n_bins = sham_data{i}.n_bins;
            if n_bins > 0
                sham_all(i, 1:min(n_bins, max_bins)) = sham_data{i}.distances(1:min(n_bins, max_bins));
            end
        end
        sham_mean = nanmean(sham_all, 1);
        sham_sem = nanstd(sham_all, 0, 1) ./ sqrt(max(sum(~isnan(sham_all), 1), 1));
        sham_n = sum(~isnan(sham_all), 1);
        
        % Ensure all arrays have correct length
        summary.SHAM_Mean = sham_mean(1:max_bins)';
        summary.SHAM_SEM = sham_sem(1:max_bins)';
        summary.SHAM_N = sham_n(1:max_bins)';
    end
    
    if ~isempty(stim_data)
        stim_all = NaN(length(stim_data), max_bins);
        for i = 1:length(stim_data)
            n_bins = stim_data{i}.n_bins;
            if n_bins > 0
                stim_all(i, 1:min(n_bins, max_bins)) = stim_data{i}.distances(1:min(n_bins, max_bins));
            end
        end
        stim_mean = nanmean(stim_all, 1);
        stim_sem = nanstd(stim_all, 0, 1) ./ sqrt(max(sum(~isnan(stim_all), 1), 1));
        stim_n = sum(~isnan(stim_all), 1);
        
        % Ensure all arrays have correct length
        summary.STIM_Mean = stim_mean(1:max_bins)';
        summary.STIM_SEM = stim_sem(1:max_bins)';
        summary.STIM_N = stim_n(1:max_bins)';
    end
    
    % Save summary
    summary_filename = sprintf('binned_locomotion_summary_%dmin.xlsx', bin_size_min);
    
    % rename *_pixels → *_mm
    summary.Properties.VariableNames = ...
        strrep(summary.Properties.VariableNames,'_pixels','_mm');

    writetable(summary, summary_filename);

    fprintf('\nSummary data saved to: %s\n', summary_filename);
end
