% plotTemporalPSDmap.m
% Draws the 2-panel time-resolved PSD map with spectrum overlay and saves .png + .fig
% Also exports raw data to Excel
%
% FEATURES:
% - Stacked bars with z-score spectrum overlay (rainbow colored segments)
% - KDE smoothing curves scaled by total peak counts
% - Standard Error of Mean (SEM) error bars across files
% - Optional Adobe .ase color palette support
% - Extended z-score range [-4σ to +4σ] with -Inf/+Inf bins to capture all peaks
% - Robust y-axis property restoration for dual-axis plots (fixes missing axis lines/labels)
%
% INPUTS:
% peakTime_STIM/SHAM    : vector of event-peak times [s]
% z_STIM/SHAM           : matching z-scores
% recDurMin             : recording length in minutes (default 30)
% binMin                : bin size in minutes (default 1)
% outDir                : folder for figures
% filenames_STIM/SHAM   : cell array of filenames corresponding to each peak
% num_stim_files/num_sham_files : number of files processed (for SEM calculation)

function plotTemporalPSDmap(peakTime_STIM, z_STIM, ...
                             peakTime_SHAM, z_SHAM, ...
                             recDurMin, binMin, outDir, ...
                             filenames_STIM, filenames_SHAM, ...
                             num_stim_files, num_sham_files)

if nargin < 4, error('Need STIM + SHAM data'); end
if nargin < 5 || isempty(recDurMin), recDurMin = 30; end
if nargin < 6 || isempty(binMin),    binMin    = 1;  end
if nargin < 7 || isempty(outDir),    outDir    = pwd; end
if nargin < 8 || isempty(filenames_STIM)
    filenames_STIM = repmat({'unknown'}, size(peakTime_STIM));
end
if nargin < 9 || isempty(filenames_SHAM)
    filenames_SHAM = repmat({'unknown'}, size(peakTime_SHAM));
end
if nargin < 10 || isempty(num_stim_files), num_stim_files = 1; end
if nargin < 11 || isempty(num_sham_files), num_sham_files = 1; end

if ~exist(outDir,'dir'), mkdir(outDir); end

% Get all user inputs in one dialog (except file inputs)
prompt = {'Enter z-score increment for spectrum overlay (e.g., 0.1, 0.25, 0.5):', ...
          'Display Mean Counts + SEM error bars? (yes/no):', ...
          'Display KDE smoothing curves? (yes/no):'};
dlgtitle = 'Analysis Parameters';
dims = [1 70];
definput = {'0.5', 'yes', 'yes'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

% Process z-score increment
if isempty(answer) || isempty(answer{1}) || isempty(str2double(answer{1})) || str2double(answer{1}) <= 0
    z_increment = 0.5; % Default to 0.5 if invalid input
    fprintf('Using default z-score increment of 0.5\n');
else
    z_increment = str2double(answer{1});
    fprintf('Using z-score increment of %g\n', z_increment);
end

% Process SEM display option
if isempty(answer) || isempty(answer{2})
    show_mean_sem = true; % Default to yes
else
    show_mean_sem = strcmpi(answer{2}, 'yes') || strcmpi(answer{2}, 'y');
end
if show_mean_sem
    fprintf('Mean counts with SEM will be displayed.\n');
else
    fprintf('Mean counts with SEM will NOT be displayed.\n');
end

% Process KDE display option
if isempty(answer) || isempty(answer{3})
    show_kde = true; % Default to yes
else
    show_kde = strcmpi(answer{3}, 'yes') || strcmpi(answer{3}, 'y');
end
if show_kde
    fprintf('KDE smoothing curves will be displayed.\n');
else
    fprintf('KDE smoothing curves will NOT be displayed.\n');
end

% Define z-score sub-bins for spectrum overlay based on user input
% Extended range to capture all z-score values including outliers
z_score_min = -4;  % Extended from -3
z_score_max = 4;   % Extended from 3
z_score_edges = [-Inf, z_score_min:z_increment:z_score_max, Inf];  % Add -Inf and Inf to capture all values

% FIX: Calculate centers correctly for the extended edge array
finite_edges = z_score_edges(2:end-1);  % Get finite edges only
finite_centers = finite_edges(1:end-1) + diff(finite_edges)/2;  % Centers of finite bins

% Create full centers array including extreme bins
z_score_centers = zeros(1, length(z_score_edges)-1);
z_score_centers(1) = z_score_min - z_increment/2;  % Center for <-4σ bin
z_score_centers(2:end-1) = finite_centers;  % Centers for finite bins
z_score_centers(end) = z_score_max + z_increment/2;  % Center for >4σ bin

num_z_bins = length(z_score_centers);

% Set up time bins
edges = 0:binMin:recDurMin;
cent  = edges(1:end-1)+binMin/2;

% Helper to bin times
toBin = @(t) discretize(t/60,edges);

% STIM - Calculate z-score distribution for each time bin
bStim = toBin(peakTime_STIM);
validStim = ~isnan(bStim);
bStim_valid = bStim(validStim);
zStim_valid = z_STIM(validStim);
filesStim_valid = filenames_STIM(validStim);
unique_stim_files = unique(filesStim_valid);
if isempty(unique_stim_files) && ~isempty(filenames_STIM) % Handle case where all peaks are outside time range but files exist
    unique_stim_files = unique(filenames_STIM);
elseif isempty(unique_stim_files) % Truly no STIM files/peaks
    unique_stim_files = {}; % Ensure it's an empty cell for consistency
end

% Initialize per-file counts for each z_sub_bin and time_bin
stimStackedDataPerFile = zeros(length(unique_stim_files), numel(edges)-1, num_z_bins);

if ~isempty(bStim_valid)
    for f_idx = 1:length(unique_stim_files)
        current_file_peaks_mask = strcmp(filesStim_valid, unique_stim_files{f_idx});
        if ~any(current_file_peaks_mask), continue; end % Skip if this file had no valid peaks

        file_bStim = bStim_valid(current_file_peaks_mask);
        file_zStim = zStim_valid(current_file_peaks_mask);

        for i = 1:numel(edges)-1 % Loop through time bins
            time_bin_mask_for_file = (file_bStim == i);
            if any(time_bin_mask_for_file)
                bin_z_scores_for_file = file_zStim(time_bin_mask_for_file);
                z_bin_indices_for_file = discretize(bin_z_scores_for_file, z_score_edges);
                for j = 1:num_z_bins % Loop through z-score sub-bins
                    stimStackedDataPerFile(f_idx, i, j) = sum(z_bin_indices_for_file == j);
                end
            end
        end
    end
    % Average across files
    stimStackedData = squeeze(mean(stimStackedDataPerFile, 1, 'omitnan'));
    % Handle potential dimension squeeze issues if num_stim_files=1 or num_z_bins=1 or num_time_bins=1
    if length(unique_stim_files) == 1 && numel(edges)-1 > 1 && num_z_bins > 1
        stimStackedData = reshape(stimStackedData, numel(edges)-1, num_z_bins);
    elseif numel(edges)-1 == 1 && num_z_bins > 1 && length(unique_stim_files) > 1
         stimStackedData = reshape(stimStackedData, 1, num_z_bins);
    elseif num_z_bins == 1 && numel(edges)-1 > 1 && length(unique_stim_files) > 1
         stimStackedData = reshape(stimStackedData, numel(edges)-1, 1);
    elseif length(unique_stim_files) == 1 && numel(edges)-1 == 1 && num_z_bins > 1
        stimStackedData = reshape(stimStackedData, 1, num_z_bins);
     elseif length(unique_stim_files) == 1 && num_z_bins == 1 && numel(edges)-1 > 1
        stimStackedData = reshape(stimStackedData, numel(edges)-1, 1);
    elseif length(unique_stim_files) > 1 && numel(edges)-1 == 1 && num_z_bins == 1
        % stimStackedData will be scalar, which is fine
    elseif length(unique_stim_files) == 1 && numel(edges)-1 == 1 && num_z_bins == 1
        % stimStackedData will be scalar, which is fine
    end
    stimStackedData(isnan(stimStackedData)) = 0; % Replace any NaN averages (e.g. no files for a bin) with 0
else
    stimStackedData = zeros(numel(edges)-1, num_z_bins);
end

% Recalculate stimCountsPerFile for SEM based on total peaks per file per time bin
stimCountsPerFile = zeros(length(unique_stim_files), numel(edges)-1);
if ~isempty(bStim_valid)
    for f_idx = 1:length(unique_stim_files)
        current_file_peaks_mask = strcmp(filesStim_valid, unique_stim_files{f_idx});
        if ~any(current_file_peaks_mask), continue; end
        file_bStim = bStim_valid(current_file_peaks_mask);
        for i = 1:numel(edges)-1 % Loop through time bins
            stimCountsPerFile(f_idx, i) = sum(file_bStim == i);
        end
    end
end

% SHAM - Calculate z-score distribution for each time bin
bSham = toBin(peakTime_SHAM);
validSham = ~isnan(bSham);
bSham_valid = bSham(validSham);
zSham_valid = z_SHAM(validSham);
filesSham_valid = filenames_SHAM(validSham);
unique_sham_files = unique(filesSham_valid);
if isempty(unique_sham_files) && ~isempty(filenames_SHAM) % Handle case where all peaks are outside time range but files exist
    unique_sham_files = unique(filenames_SHAM);
elseif isempty(unique_sham_files) % Truly no SHAM files/peaks
    unique_sham_files = {}; % Ensure it's an empty cell for consistency
end

% Initialize per-file counts for each z_sub_bin and time_bin
shamStackedDataPerFile = zeros(length(unique_sham_files), numel(edges)-1, num_z_bins);

if ~isempty(bSham_valid)
    for f_idx = 1:length(unique_sham_files)
        current_file_peaks_mask = strcmp(filesSham_valid, unique_sham_files{f_idx});
        if ~any(current_file_peaks_mask), continue; end % Skip if this file had no valid peaks

        file_bSham = bSham_valid(current_file_peaks_mask);
        file_zSham = zSham_valid(current_file_peaks_mask);

        for i = 1:numel(edges)-1 % Loop through time bins
            time_bin_mask_for_file = (file_bSham == i);
            if any(time_bin_mask_for_file)
                bin_z_scores_for_file = file_zSham(time_bin_mask_for_file);
                z_bin_indices_for_file = discretize(bin_z_scores_for_file, z_score_edges);
                for j = 1:num_z_bins % Loop through z-score sub-bins
                    shamStackedDataPerFile(f_idx, i, j) = sum(z_bin_indices_for_file == j);
                end
            end
        end
    end
    % Average across files
    shamStackedData = squeeze(mean(shamStackedDataPerFile, 1, 'omitnan'));
    % Handle potential dimension squeeze issues if num_sham_files=1 or num_z_bins=1 or num_time_bins=1
    if length(unique_sham_files) == 1 && numel(edges)-1 > 1 && num_z_bins > 1
        shamStackedData = reshape(shamStackedData, numel(edges)-1, num_z_bins);
    elseif numel(edges)-1 == 1 && num_z_bins > 1 && length(unique_sham_files) > 1
         shamStackedData = reshape(shamStackedData, 1, num_z_bins);
    elseif num_z_bins == 1 && numel(edges)-1 > 1 && length(unique_sham_files) > 1
         shamStackedData = reshape(shamStackedData, numel(edges)-1, 1);
    elseif length(unique_sham_files) == 1 && numel(edges)-1 == 1 && num_z_bins > 1
        shamStackedData = reshape(shamStackedData, 1, num_z_bins);
    elseif length(unique_sham_files) == 1 && num_z_bins == 1 && numel(edges)-1 > 1
        shamStackedData = reshape(shamStackedData, numel(edges)-1, 1);
    elseif length(unique_sham_files) > 1 && numel(edges)-1 == 1 && num_z_bins == 1
        % shamStackedData will be scalar, which is fine
    elseif length(unique_sham_files) == 1 && numel(edges)-1 == 1 && num_z_bins == 1
        % shamStackedData will be scalar, which is fine
    end
    shamStackedData(isnan(shamStackedData)) = 0; % Replace any NaN averages (e.g. no files for a bin) with 0
else
    shamStackedData = zeros(numel(edges)-1, num_z_bins);
end

% Recalculate shamCountsPerFile for SEM based on total peaks per file per time bin
shamCountsPerFile = zeros(length(unique_sham_files), numel(edges)-1);
if ~isempty(bSham_valid)
    for f_idx = 1:length(unique_sham_files)
        current_file_peaks_mask = strcmp(filesSham_valid, unique_sham_files{f_idx});
        if ~any(current_file_peaks_mask), continue; end
        file_bSham = bSham_valid(current_file_peaks_mask);
        for i = 1:numel(edges)-1 % Loop through time bins
            shamCountsPerFile(f_idx, i) = sum(file_bSham == i);
        end
    end
end

% Calculate SEM for STIM and SHAM
if size(stimCountsPerFile, 1) > 1
    stimSEM = std(stimCountsPerFile, 0, 1, 'omitnan') ./ sqrt(size(stimCountsPerFile, 1));
    stimMeanCounts = mean(stimCountsPerFile, 1, 'omitnan');
else
    stimSEM = zeros(1, size(stimCountsPerFile, 2));
    stimMeanCounts = stimCountsPerFile; % If only one file, mean is the file itself, SEM is 0
end

if size(shamCountsPerFile, 1) > 1
    shamSEM = std(shamCountsPerFile, 0, 1, 'omitnan') ./ sqrt(size(shamCountsPerFile, 1));
    shamMeanCounts = mean(shamCountsPerFile, 1, 'omitnan');
else
    shamSEM = zeros(1, size(shamCountsPerFile, 2));
    shamMeanCounts = shamCountsPerFile; % If only one file, mean is the file itself, SEM is 0
end

% Export raw data to Excel
% Create data table for Excel export
stimData = table();
stimData.Filename = filenames_STIM(:);
stimData.Time_min = peakTime_STIM(:)/60;  % Convert to minutes
stimData.ZScore = z_STIM(:);
stimData.Condition = repmat({'STIM'}, length(peakTime_STIM), 1);

shamData = table();
shamData.Filename = filenames_SHAM(:);
shamData.Time_min = peakTime_SHAM(:)/60;  % Convert to minutes
shamData.ZScore = z_SHAM(:);
shamData.Condition = repmat({'SHAM'}, length(peakTime_SHAM), 1);

% Combine all data
allData = [stimData; shamData];

% Save to Excel
excelFile = fullfile(outDir, 'Temporal_PSD_RawData.xlsx');
writetable(allData, excelFile);
fprintf('Raw data exported to: %s\n', excelFile);

% Calculate KDE smoothing for probability density
% Create smooth time range for KDE evaluation
kde_time_range = linspace(0, recDurMin, 200);  % 200 points for smooth curves

% Calculate KDE for STIM with scaled output
if ~isempty(peakTime_STIM)
    total_stim_peaks = length(peakTime_STIM);
    try
        [stim_kde_density, ~] = ksdensity(peakTime_STIM/60, kde_time_range, 'Bandwidth', binMin);
        % Scale KDE by total peak count to get actual peak density
        stimKDESmooth = stim_kde_density * total_stim_peaks;
    catch
        % Fallback if KDE fails
        stimKDESmooth = zeros(size(kde_time_range));
    end
else
    stimKDESmooth = zeros(size(kde_time_range));
end

% Calculate KDE for SHAM with scaled output
if ~isempty(peakTime_SHAM)
    total_sham_peaks = length(peakTime_SHAM);
    try
        [sham_kde_density, ~] = ksdensity(peakTime_SHAM/60, kde_time_range, 'Bandwidth', binMin);
        % Scale KDE by total peak count to get actual peak density
        shamKDESmooth = sham_kde_density * total_sham_peaks;
    catch
        % Fallback if KDE fails
        shamKDESmooth = zeros(size(kde_time_range));
    end
else
    shamKDESmooth = zeros(size(kde_time_range));
end

% Plotting
try
    cmap_levels = 256;
    cmap = jet(cmap_levels);   % Default fallback, not used for stacked bars
catch
    cmap_levels = 256;
    cmap = parula(cmap_levels); % Default fallback
end

figure('Position',[100 100 1400 500],'Name','Temporal PSD map');

% Ask user if they want to use custom Adobe .ase color palette
use_custom_palette = questdlg('Do you want to use a custom Adobe .ase color palette for z-scores?', ...
                              'Color Palette Selection', ...
                              'Yes', 'No (use default jet)', 'No (use default jet)');

if strcmp(use_custom_palette, 'Yes')
    [ase_filename, ase_pathname] = uigetfile('*.ase', 'Select Adobe .ase color palette file');
    if ~isequal(ase_filename, 0)
        try
            ase_colors = readASEFile(fullfile(ase_pathname, ase_filename));
            if ~isempty(ase_colors)
                % MODIFICATION 2: Flip/inverse the ASE color palette
                ase_colors = flipud(ase_colors);  % Flip the color order
                
                if size(ase_colors, 1) >= num_z_bins
                    indices = round(linspace(1, size(ase_colors, 1), num_z_bins));
                    z_cmap = ase_colors(indices, :);
                else
                    z_cmap = interpolateColors(ase_colors, num_z_bins);
                end
                fprintf('Successfully loaded custom color palette with %d colors for %d z-bins (INVERTED)\n', size(ase_colors, 1), num_z_bins);
            else
                fprintf('Warning: Could not read colors from .ase file, using default jet colormap\n');
                z_cmap = jet(num_z_bins);
            end
        catch err
            fprintf('Warning: Error reading .ase file (%s), using default jet colormap\n', err.message);
            z_cmap = jet(num_z_bins);
        end
    else
        fprintf('No .ase file selected, using default jet colormap\n');
        z_cmap = jet(num_z_bins);
    end
else
    z_cmap = jet(num_z_bins);
end

% Calculate maxCount for y-axis scaling
if show_mean_sem
    max_stim_val = max(stimMeanCounts + stimSEM);
    max_sham_val = max(shamMeanCounts + shamSEM);
    % Ensure stimMeanCounts and shamMeanCounts are row vectors for sum if they are single values
    if isscalar(stimMeanCounts), stimMeanCounts_for_sum = stimMeanCounts; else, stimMeanCounts_for_sum = sum(stimMeanCounts,2); end;
    if isscalar(shamMeanCounts), shamMeanCounts_for_sum = shamMeanCounts; else, shamMeanCounts_for_sum = sum(shamMeanCounts,2); end;
    maxCount = max([max_stim_val, max_sham_val, max(sum(stimStackedData,2)), max(sum(shamStackedData,2))]);
else
    maxCount = max([max(sum(stimStackedData,2)), max(sum(shamStackedData,2))]);
end
if isempty(maxCount) || isnan(maxCount) || maxCount == 0, maxCount = 1; end % Prevent zero or NaN ylim.

% MODIFICATION 3: Determine y-axis tick interval based on maxCount
if maxCount > 20
    y_tick_interval = 20;  % Use 20 when max exceeds 20
else
    y_tick_interval = 10;  % Use 10 for smaller values
end

% Calculate maxKDESmooth for KDE y-axis scaling
maxKDESmooth = max([max(stimKDESmooth), max(shamKDESmooth)]);
if isempty(maxKDESmooth) || isnan(maxKDESmooth) || maxKDESmooth == 0, maxKDESmooth = 1;
end

% STIM panel
subplot(1,2,2);
yyaxis left  % Left y-axis for peak counts
hold on; box on; grid on;

% Plot stacked bars for STIM with spectrum overlay
bh_stim = bar(cent, stimStackedData, 'stacked', 'EdgeColor', 'k', 'LineWidth', 0.5);
% Set colors for each z-score range
for j = 1:num_z_bins
    bh_stim(j).FaceColor = z_cmap(j, :);
end

% Add error bars for SEM if requested
if show_mean_sem
    hold on;
    errorbar(cent, stimMeanCounts, stimSEM, 'k.', 'LineWidth', 1.5, 'CapSize', 3);
    fprintf('STIM: SEM error bars plotted.\n'); % Debugging line
end

ylabel('Average number of peaks / file');
ylim([0, maxCount * 1.1]);
% MODIFICATION 3: Set y-axis ticks based on determined interval
yticks(0:y_tick_interval:ceil(maxCount * 1.1));
ax1 = gca;
ax1.YColor = 'k';

% EXPLICIT Y-AXIS FIX IMPLEMENTATION:
% Purpose: Ensures left Y-axis properties (line, labels, ticks) are consistently correct
% Cause: Switching axes using yyaxis left/right can overwrite properties unintentionally
% Solution: Separate axis handles (ax_left/ax_right) and explicitly reset properties

if show_kde || show_mean_sem
    % Normal plotting with KDE or SEM enabled
    yyaxis left
    ax1_left = gca;
    ax1_left.YColor = 'k';
    ylabel('Average number of peaks / file');
    ylim([0, maxCount * 1.1]);
    yticks(0:y_tick_interval:ceil(maxCount * 1.1));
    
    if show_kde
        yyaxis right  % Right y-axis for KDE smoothing
        ax1_right = gca;
        plot(kde_time_range, stimKDESmooth, 'r-', 'LineWidth', 2);
        ylabel('Peak density (KDE scaled)');
        ylim([0, maxKDESmooth * 1.1]);
        ax1_right.YColor = 'r';
        % Ensure box is on for right axis
        ax1_right.Box = 'on';
    else
        % SEM only, hide right axis but keep left axis properties
        yyaxis right
        ax1_right = gca;
        ax1_right.YColor = 'none';
        ax1_right.YTick = [];
        ax1_right.YTickLabel = [];
        ylabel('');
        % Ensure box is on even when right axis is hidden
        ax1_right.Box = 'on';
        yyaxis left  % Return to left axis
    end
else
    % ----- Robust Comprehensive Y-axis Display Error Fix -----
    % Scenario: Both KDE curves and SEM bars are disabled
    yyaxis right
    ax_right = gca;
    ax_right.YColor = 'none';
    ax_right.YTick = [];
    ax_right.YTickLabel = [];
    ylabel('');
    % Ensure box is on even when right axis is hidden
    ax_right.Box = 'on';

    yyaxis left
    ax_left = gca;
    ax_left.YColor = 'k';
    ylabel('Average number of peaks / file', 'Color', 'k');

    maxY = ceil(maxCount * 1.1);
    ylim([0, maxY]);
    yticks(0:y_tick_interval:maxY);

    xlim([0 recDurMin]); % Explicit x-axis synchronization

    ax_left.XColor = 'k';
    ax_left.Box = 'on';
    ax_left.GridLineStyle = '-';
    grid on;
    hold on;
end

yyaxis left  % Switch back to left for other properties
title('STIM'); xlabel('Time [min]');
xlim([0 recDurMin]);

% FIX 1: Ensure complete box is drawn
set(gca, 'Box', 'on');
% Force redraw of axes to ensure all borders are visible
drawnow;

% SHAM panel
subplot(1,2,1);
yyaxis left  % Left y-axis for peak counts
hold on; box on; grid on;

% Plot stacked bars for SHAM with spectrum overlay
bh_sham = bar(cent, shamStackedData, 'stacked', 'EdgeColor', 'k', 'LineWidth', 0.5);
% Set colors for each z-score range
for j = 1:num_z_bins
    bh_sham(j).FaceColor = z_cmap(j, :);
end

% Add error bars for SEM if requested
if show_mean_sem
    hold on;
    errorbar(cent, shamMeanCounts, shamSEM, 'k.', 'LineWidth', 1.5, 'CapSize', 3);
    fprintf('SHAM: SEM error bars plotted.\n'); % Debugging line
end

ylabel('Average number of peaks / file');
ylim([0, maxCount * 1.1]);
% MODIFICATION 3: Set y-axis ticks based on determined interval
yticks(0:y_tick_interval:ceil(maxCount * 1.1));
ax2 = gca;
ax2.YColor = 'k';

% EXPLICIT Y-AXIS FIX IMPLEMENTATION (SHAM Panel):
% Purpose: Ensures left Y-axis properties (line, labels, ticks) are consistently correct
% Cause: Switching axes using yyaxis left/right can overwrite properties unintentionally
% Solution: Separate axis handles (ax_left/ax_right) and explicitly reset properties

if show_kde || show_mean_sem
    % Normal plotting with KDE or SEM enabled
    yyaxis left
    ax2_left = gca;
    ax2_left.YColor = 'k';
    ylabel('Average number of peaks / file');
    ylim([0, maxCount * 1.1]);
    yticks(0:y_tick_interval:ceil(maxCount * 1.1));
    
    if show_kde
        yyaxis right  % Right y-axis for KDE smoothing
        ax2_right = gca;
        plot(kde_time_range, shamKDESmooth, 'r-', 'LineWidth', 2);
        ylabel('Peak density (KDE scaled)');
        ylim([0, maxKDESmooth * 1.1]);
        ax2_right.YColor = 'r';
        % Ensure box is on for right axis
        ax2_right.Box = 'on';
    else
        % SEM only, hide right axis but keep left axis properties
        yyaxis right
        ax2_right = gca;
        ax2_right.YColor = 'none';
        ax2_right.YTick = [];
        ax2_right.YTickLabel = [];
        ylabel('');
        % Ensure box is on even when right axis is hidden
        ax2_right.Box = 'on';
        yyaxis left  % Return to left axis
    end
else
    % ----- Robust Comprehensive Y-axis Display Error Fix -----
    % Scenario: Both KDE curves and SEM bars are disabled
    yyaxis right
    ax_right = gca;
    ax_right.YColor = 'none';
    ax_right.YTick = [];
    ax_right.YTickLabel = [];
    ylabel('');
    % Ensure box is on even when right axis is hidden
    ax_right.Box = 'on';

    yyaxis left
    ax_left = gca;
    ax_left.YColor = 'k';
    ylabel('Average number of peaks / file', 'Color', 'k');

    maxY = ceil(maxCount * 1.1);
    ylim([0, maxY]);
    yticks(0:y_tick_interval:maxY);

    xlim([0 recDurMin]); % Explicit x-axis synchronization

    ax_left.XColor = 'k';
    ax_left.Box = 'on';
    ax_left.GridLineStyle = '-';
    grid on;
    hold on;
end

yyaxis left  % Switch back to left for other properties
title('SHAM'); xlabel('Time [min]');
xlim([0 recDurMin]);

% FIX 1: Ensure complete box is drawn
set(gca, 'Box', 'on');
% Force redraw of axes to ensure all borders are visible
drawnow;

% Shared colorbar
cb = colorbar('Position',[0.92 0.13 0.02 0.74]);
colormap(z_cmap);
% Use finite z-score range for caxis (exclude -Inf and Inf)
finite_edges = z_score_edges(~isinf(z_score_edges));
caxis([min(finite_edges) max(finite_edges)]);

% MODIFICATION 1: Set colorbar ticks to every 2.0 sigma for the finite range
tick_values = z_score_min:2:z_score_max;  % Changed from 1 to 2
% Ensure ticks are within the actual caxis limits and are unique
tick_values = unique(max(min(finite_edges), min(max(finite_edges), tick_values)));
cb.Ticks = tick_values;

ylabel(cb,'z-score range (\sigma)','FontWeight','bold');
sgtitle(sprintf('Temporal distribution of detected peaks (bin = %d min)', binMin));

% Save figures
saveas(gcf,fullfile(outDir,'Temporal_PSDmap.png'));
saveas(gcf,fullfile(outDir,'Temporal_PSDmap.fig'));

% Print summary statistics
fprintf('\nTemporal PSD map generated:\n');
if show_mean_sem
    fprintf('STIM: %d total peaks across %d minutes (%d files, Mean+SEM displayed)\n', sum(sum(stimStackedData)), recDurMin, length(unique_stim_files));
    fprintf('SHAM: %d total peaks across %d minutes (%d files, Mean+SEM displayed)\n', sum(sum(shamStackedData)), recDurMin, length(unique_sham_files));
else
    fprintf('STIM: %d total peaks across %d minutes (%d files, Mean+SEM NOT displayed)\n', sum(sum(stimStackedData)), recDurMin, length(unique_stim_files));
    fprintf('SHAM: %d total peaks across %d minutes (%d files, Mean+SEM NOT displayed)\n', sum(sum(shamStackedData)), recDurMin, length(unique_sham_files));
end
if show_kde
    fprintf('Smoothing: KDE with bandwidth = %.1f min, scaled by peak counts\n', binMin);
else
    fprintf('Smoothing: KDE curves disabled by user\n');
end

end

%% Helper Functions

function colors = readASEFile(filename)
    % Read Adobe Swatch Exchange (.ase) file
    % Returns an Nx3 matrix of RGB colors (values 0-1)
    
    colors = [];
    
    try
        fid = fopen(filename, 'rb');
        if fid == -1
            error('Could not open file');
        end
        
        % Read ASE header
        signature = fread(fid, 4, '*char')';
        if ~strcmp(signature, 'ASEF')
            fclose(fid);
            error('Not a valid ASE file');
        end
        
        version_major = fread(fid, 1, 'uint16', 0, 'ieee-be');
        version_minor = fread(fid, 1, 'uint16', 0, 'ieee-be');
        num_blocks = fread(fid, 1, 'uint32', 0, 'ieee-be');
        
        colors = [];
        
        for i = 1:num_blocks
            % Read block header
            block_type = fread(fid, 1, 'uint16', 0, 'ieee-be');
            block_length = fread(fid, 1, 'uint32', 0, 'ieee-be');
            
            if block_type == 1  % Color entry
                % Read name length and name
                name_length = fread(fid, 1, 'uint16', 0, 'ieee-be');
                if name_length > 0
                    name = fread(fid, name_length*2, '*char')';  % UTF-16, skip for now
                end
                
                % Read color model
                color_model = fread(fid, 4, '*char')';
                
                if strcmp(color_model, 'RGB ')
                    % Read RGB values
                    r = fread(fid, 1, 'single', 0, 'ieee-be');
                    g = fread(fid, 1, 'single', 0, 'ieee-be');
                    b = fread(fid, 1, 'single', 0, 'ieee-be');
                    
                    % Add to colors array
                    colors = [colors; r, g, b];
                else
                    % Skip non-RGB colors
                    remaining_bytes = block_length - 6 - name_length*2 - 4;
                    if remaining_bytes > 0
                        fread(fid, remaining_bytes, 'uint8');
                    end
                end
                
                % Read color type (spot/process)
                color_type = fread(fid, 1, 'uint16', 0, 'ieee-be');
                
            else
                % Skip non-color blocks
                if block_length > 0
                    fread(fid, block_length, 'uint8');
                end
            end
        end
        
        fclose(fid);
        
        % Ensure colors are in valid range [0,1]
        colors = max(0, min(1, colors));
        
    catch err
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        warning('Error reading ASE file: %s', err.message);
        colors = [];
    end
end

function interpolated_colors = interpolateColors(colors, num_colors)
    % Interpolate colors to get the desired number of colors
    % colors: Nx3 matrix of RGB values
    % num_colors: desired number of output colors
    
    if isempty(colors)
        interpolated_colors = jet(num_colors);
        return;
    end
    
    if size(colors, 1) == 1
        % Only one color, replicate it
        interpolated_colors = repmat(colors, num_colors, 1);
        return;
    end
    
    % Create interpolation indices
    original_indices = 1:size(colors, 1);
    target_indices = linspace(1, size(colors, 1), num_colors);
    
    % Interpolate each color channel
    interpolated_colors = zeros(num_colors, 3);
    for ch = 1:3
        interpolated_colors(:, ch) = interp1(original_indices, colors(:, ch), target_indices, 'linear');
    end
    
    % Ensure values are in valid range
    interpolated_colors = max(0, min(1, interpolated_colors));
end
