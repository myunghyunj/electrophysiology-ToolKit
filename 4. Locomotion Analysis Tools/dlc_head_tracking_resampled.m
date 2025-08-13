% dlc_head_tracking_resampled.m
% Track head distance from DLC exports with resampling from 30fps to 3fps
% Calculates Euclidean distances and timestamps the data

%% USER PARAMETERS
likelihood_thresh = 0.90;       % confidence cutoff (same as batch_DLC_locomotion.m)
fps_original = 30;              % original video fps
fps_target = 3;                 % target fps after resampling
resample_factor = fps_original / fps_target;  % downsample by factor of 10

fprintf('\n=== DLC Head Tracking with Resampling ===\n');
fprintf('Original FPS: %d | Target FPS: %d | Downsample factor: %d\n', ...
    fps_original, fps_target, resample_factor);
fprintf('Processing CSV files from raw folder...\n\n');

%% SELECT CSV FILE OR PROCESS ALL
choice = input('Process [1] single file or [2] all CSV files? Enter 1 or 2: ');

if choice == 1
    % Single file processing
    [filename, filepath] = uigetfile('raw/*.csv', 'Select DLC CSV file');
    if filename == 0
        error('No file selected');
    end
    files = dir(fullfile(filepath, filename));
else
    % Batch processing
    files = dir('raw/*.csv');
    if isempty(files)
        error('No CSV files found in raw folder');
    end
end

%% PROCESS EACH FILE
for f = 1:numel(files)
    fname = files(f).name;
    fpath = fullfile('raw', fname);
    
    fprintf('Processing: %s\n', fname);
    
    % Read CSV file - skip first 3 header rows
    data = readmatrix(fpath, 'NumHeaderLines', 3);
    
    % Read header to get column structure
    fid = fopen(fpath, 'r');
    line1 = fgetl(fid);  % scorer
    line2 = fgetl(fid);  % bodyparts
    line3 = fgetl(fid);  % coords
    fclose(fid);
    
    % Parse bodyparts to find head columns
    bodyparts = strsplit(line2, ',');
    coords = strsplit(line3, ',');
    
    % Find head x, y, likelihood columns (1-indexed, accounting for frame column)
    head_idx = find(strcmp(bodyparts, 'head'));
    if isempty(head_idx)
        warning('No head marker found in %s', fname);
        continue;
    end
    
    % First occurrence of 'head' corresponds to x coordinate
    head_x_col = head_idx(1);
    head_y_col = head_x_col + 1;
    head_like_col = head_x_col + 2;
    
    % Extract data
    frames = data(:, 1);
    head_x = data(:, head_x_col);
    head_y = data(:, head_y_col);
    head_likelihood = data(:, head_like_col);
    
    % Apply likelihood threshold
    valid_idx = head_likelihood >= likelihood_thresh;
    head_x_filtered = head_x;
    head_y_filtered = head_y;
    head_x_filtered(~valid_idx) = NaN;
    head_y_filtered(~valid_idx) = NaN;
    
    % Fill missing values using pchip interpolation (same as batch_DLC_locomotion.m)
    head_x_interp = fillmissing(head_x_filtered, 'pchip');
    head_y_interp = fillmissing(head_y_filtered, 'pchip');
    
    % Resample data from 30fps to 3fps (take every 10th frame)
    resample_indices = 1:resample_factor:length(head_x_interp);
    head_x_resampled = head_x_interp(resample_indices);
    head_y_resampled = head_y_interp(resample_indices);
    frames_resampled = frames(resample_indices);
    
    % Calculate Euclidean distances between consecutive resampled frames
    distances = sqrt(diff(head_x_resampled).^2 + diff(head_y_resampled).^2);
    
    % Create timestamps (in seconds) for resampled data
    time_seconds = (frames_resampled - frames_resampled(1)) / fps_original;
    
    % Create results table
    % Note: distances array is one element shorter than positions
    results = table();
    results.Frame = frames_resampled(1:end-1);
    results.Time_sec = time_seconds(1:end-1);
    results.Head_X = head_x_resampled(1:end-1);
    results.Head_Y = head_y_resampled(1:end-1);
    results.Distance_pixels = distances;
    results.Speed_pixels_per_sec = distances * fps_target;  % Convert to pixels/second
    
    % Calculate summary statistics
    total_distance = nansum(distances);
    total_time = time_seconds(end);
    avg_speed = total_distance / total_time;
    max_speed = max(distances) * fps_target;
    
    % Save results to Excel file
    [~, basename, ~] = fileparts(fname);
    output_fname = sprintf('%s_head_tracking_3fps.xlsx', basename);
    writetable(results, output_fname);
    
    % Display summary
    fprintf('  Original frames: %d | Resampled frames: %d\n', ...
        length(frames), length(frames_resampled));
    fprintf('  Total distance: %.2f pixels\n', total_distance);
    fprintf('  Total time: %.2f seconds (%.2f minutes)\n', total_time, total_time/60);
    fprintf('  Average speed: %.2f pixels/sec\n', avg_speed);
    fprintf('  Maximum speed: %.2f pixels/sec\n', max_speed);
    fprintf('  Output saved to: %s\n\n', output_fname);
    
    % Create trajectory plot for all files
    fig = figure('Name', sprintf('Head Trajectory - %s', basename));
    subplot(2,1,1);
    
    % Use different colors for SHAM vs STIM
    if contains(upper(fname), 'SHAM')
        plot_color = 'b-';  % Blue for SHAM
    else
        plot_color = 'r-';  % Red for STIM
    end
    
    plot(head_x_resampled, head_y_resampled, plot_color, 'LineWidth', 1.5);
    hold on;
    plot(head_x_resampled(1), head_y_resampled(1), 'go', 'MarkerSize', 8, 'LineWidth', 2);
    plot(head_x_resampled(end), head_y_resampled(end), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
    xlabel('X (pixels)'); ylabel('Y (pixels)');
    title(sprintf('Head Trajectory (3fps resampled) - %s', basename), 'Interpreter', 'none');
    legend('Trajectory', 'Start', 'End', 'Location', 'best');
    axis equal; grid on;
    
    subplot(2,1,2);
    plot(time_seconds(1:end-1), results.Speed_pixels_per_sec, 'k-', 'LineWidth', 1);
    xlabel('Time (seconds)'); ylabel('Speed (pixels/sec)');
    title('Head Speed Over Time');
    grid on;
    
    % Save figure
    saveas(fig, sprintf('%s_head_analysis_3fps.png', basename));
    
    % Close figure to prevent memory issues with many files
    if choice == 2  % Only close in batch mode
        close(fig);
    end
end

fprintf('=== Processing Complete ===\n');
