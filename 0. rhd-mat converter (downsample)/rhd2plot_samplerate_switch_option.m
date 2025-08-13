% convert_rhd_to_mat (downsamples if over sampled)
%   Reads an Intan *.rhd file, extracts amplifier_data → raw_data,
%   prints the original sampling rate, down‑samples to 1 kHz when needed,
%   and saves raw_data plus Fs_final into a *.mat with the same base‑name.
%
%   Requirements: read_Intan_RHD2000_file.m must be on the MATLAB path.
%                 This script should be renamed (e.g., to rhd2plot_samplerate.m) to avoid issues with spaces in the filename.

%% 1) Clear leftovers from previous runs
vars_to_clear = {'filename','path','notes','frequency_parameters', ...
                 'reference_channel','amplifier_channels','amplifier_data', ...
                 't_amplifier','spike_triggers','aux_input_channels', ...
                 'aux_input_data','t_aux_input','supply_voltage_channels', ...
                 'supply_voltage_data','t_supply_voltage','board_adc_channels', ...
                 'board_adc_data','t_board_adc','board_dig_in_channels', ...
                 'board_dig_in_data','t_dig','board_dig_out_channels', ...
                 'board_dig_out_data','temp_sensor_data','t_temp_sensor'};
evalin('base',['clearvars ', strjoin(vars_to_clear,' ')]);

%% 2) Read *.rhd
disp('Select an RHD file …');
evalin('base','read_Intan_RHD2000_file');

if ~(evalin('base','exist(''filename'',''var'')') && evalin('base','~isempty(filename)'))
    disp('No file picked. Exit.'); return; end

fname = evalin('base','filename');
fpath = evalin('base','path');
fprintf('Processing: %s\n', fullfile(fpath,fname));

if ~evalin('base','exist(''frequency_parameters'',''var'')')
    disp('frequency_parameters missing. Exit.'); return; end
freq = evalin('base','frequency_parameters');
Fs_orig = freq.amplifier_sample_rate;

if ~evalin('base','exist(''amplifier_data'',''var'')')
    fprintf('amplifier_data missing. Fs = %.2f Hz. Exit.\n', Fs_orig); return; end
raw_data = evalin('base','amplifier_data');
amplifier_channels = evalin('base','amplifier_channels');

%% 3) Optional channel‑row reorder
answ = input('Reorder channel rows before saving? Y/N [N]: ','s');
if strcmpi(answ,'Y')
    fprintf('Current order:\n');
    for k = 1:numel(amplifier_channels)
        fprintf('%3d: %s\n',k,amplifier_channels(k).native_channel_name);
    end
    newOrder = input(sprintf('Enter new order as index vector [1‑%d]: ',numel(amplifier_channels)));
    if numel(newOrder) ~= size(raw_data,1) || ~all(sort(newOrder)==1:size(raw_data,1))
        error('Invalid vector. Abort.');
    end
    raw_data = raw_data(newOrder,:);
    amplifier_channels = amplifier_channels(newOrder);
    fprintf('Rows reordered.\n');
end

%% 4) Down‑sample when needed
fprintf('Original Fs: %.2f Hz\n',Fs_orig);
if Fs_orig > 1000
    q = round(Fs_orig/1000); if q==0, q=1; end
    fprintf('Down‑sampling ×%d to %.2f Hz …\n',q,Fs_orig/q);
    raw_data = downsample(raw_data.',q).';
    Fs_final = Fs_orig/q;
else
    Fs_final = Fs_orig;
end

%% 5) Save
[~,base] = fileparts(fname);
outfile = fullfile(fpath,[base '.mat']);
fprintf('Saving → %s …\n',outfile);
save(outfile,'raw_data','Fs_final','-v7.3');

fprintf('Done! %d ch × %d samples @ %.2f Hz\n',size(raw_data,1),size(raw_data,2),Fs_final);
