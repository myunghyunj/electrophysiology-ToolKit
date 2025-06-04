% convert_rhd_to_mat
%   Reads an Intan *.rhd file, extracts amplifier_data → raw_data,
%   prints the original sampling rate, down‑samples to 1 kHz when needed,
%   and saves raw_data plus Fs_final into a *.mat with the same base‑name.
%
%   Requirements: read_Intan_RHD2000_file.m must be on the MATLAB path.
%                 This script should be renamed (e.g., to rhd2plot_samplerate.m) to avoid issues with spaces in the filename.

    %% 1) Clear potentially leftover variables from previous runs in base workspace
    % These are variables that read_Intan_RHD2000_file.m is known to create.
    vars_to_clear = {'filename', 'path', 'notes', 'frequency_parameters', ...
                     'reference_channel', 'amplifier_channels', 'amplifier_data', ...
                     't_amplifier', 'spike_triggers', 'aux_input_channels', ...
                     'aux_input_data', 't_aux_input', 'supply_voltage_channels', ...
                     'supply_voltage_data', 't_supply_voltage', 'board_adc_channels', ...
                     'board_adc_data', 't_board_adc', 'board_dig_in_channels', ...
                     'board_dig_in_data', 't_dig', 'board_dig_out_channels', ...
                     'board_dig_out_data', 'temp_sensor_data', 't_temp_sensor'};
    evalin('base', ['clearvars ', strjoin(vars_to_clear, ' ')]);

    %% 2) Read the RHD file using read_Intan_RHD2000_file
    % This script will prompt for file selection and populate variables in the base workspace.
    disp('Please select an RHD file in the dialog opened by read_Intan_RHD2000_file.m...');
    evalin('base','read_Intan_RHD2000_file');

    %% 3) Fetch data, sampling rate, and actual file info from base workspace
    % Check if read_Intan_RHD2000_file successfully selected a file and set variables
    if ~(evalin('base','exist(''filename'',''var'')') && evalin('base','~isempty(filename)'))
        disp('File selection was cancelled or failed within read_Intan_RHD2000_file. Exiting.');
        return;
    end

    actual_rhd_file_name = evalin('base','filename');
    actual_rhd_file_path = evalin('base','path');
    fprintf('Processing RHD file: %s\n', fullfile(actual_rhd_file_path, actual_rhd_file_name));

    if ~evalin('base','exist(''frequency_parameters'',''var'')')
        disp('Error: frequency_parameters not found after running read_Intan_RHD2000_file. Exiting.');
        return;
    end
    freq_param  = evalin('base','frequency_parameters');
    Fs_original = freq_param.amplifier_sample_rate;

    % amplifier_data is crucial for processing.
    if ~evalin('base','exist(''amplifier_data'',''var'')')
        fprintf('Warning: amplifier_data not found. The file might be a header-only file or lack amplifier data.\n');
        fprintf('Original sampling rate from header: %.2f Hz. Cannot proceed with data processing.\n', Fs_original);
        return;
    end
    raw_data    = evalin('base','amplifier_data');
    
    %% 4) Display the native sampling rate
    fprintf('Original sampling rate: %.2f Hz\n', Fs_original);

    %% 5) Down‑sample if necessary
    if Fs_original > 1000
        q = round(Fs_original/1000);                 % integer factor
        % Ensure q is at least 1 to prevent errors with downsample if Fs_original is, e.g., 1000.1 Hz
        if q == 0; q = 1; end
        fprintf('Down-sampling by %d times to %.2f Hz...\n', q, Fs_original/q);
        raw_data  = downsample(raw_data.', q).';     % operate column‑wise
        Fs_final  = Fs_original/q;
    else
        Fs_final  = Fs_original;
    end

    %% 6) Save
    [~, base_name] = fileparts(actual_rhd_file_name); % Use the actual file name read
    outfile   = fullfile(actual_rhd_file_path, [base_name '.mat']); % Use the actual path from the reader
    fprintf('Saving to "%s"...\n', outfile);
    save(outfile,'raw_data','Fs_final','-v7.3');

    fprintf('Done!  %d channels × %d samples @ %.2f Hz\n', ...
            size(raw_data,1), size(raw_data,2), Fs_final);