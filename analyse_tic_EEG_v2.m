function analyse_tic_EEG_v2(dataFile)
% ANALYSE_TIC_EEG_V2  Complete tic EEG analysis for Intan recordings.
%   ANALYSE_TIC_EEG_V2(DATAFILE) processes epidural screw EEG recorded with
%   six Intan inputs (5 active + 1 ground). DATAFILE can be a *.rhd file or
%   a MATLAB *.mat file containing variables 'amplifier_data' and 'fs'.
%   Results are saved to 'results.mat' and a summary figure is written as
%   'summary_<date>.png'.
%
%   The output structure includes pre‑processed EEG, burst timings, power
%   spectra, coherence, and phase metrics.
%
%   Example:
%       analyse_tic_EEG_v2('session1.rhd');
%
%   Requirements: Signal Processing Toolbox.

if nargin < 1 || isempty(dataFile)
    [f,p] = uigetfile({'*.rhd;*.mat','RHD or MAT'}, 'Select data file');
    if isequal(f,0)
        error('No file selected.');
    end
    dataFile = fullfile(p,f);
end

[eegRaw, fs] = load_intan(dataFile);

[eegProc, t, fs] = preprocess(eegRaw, fs);

bursts = detectBursts(eegProc(:,1:4), fs);

[coherence, phaseLag] = computeCoherence(eegProc, bursts, fs);

bandPower = computeBandPower(eegProc, fs);

plotSummary(t, eegProc, bursts, coherence, fs);

save('results.mat','fs','t','eegProc','bursts','bandPower','coherence','phaseLag');
end

%% Local Functions
function [data, fs] = load_intan(fileName)
%LOAD_INTAN Load Intan RHD or MAT data.
%   [DATA, FS] = LOAD_INTAN(FILENAME) reads amplifier data and sampling
%   rate from either a .mat file or a .rhd file. DATA is returned as
%   samples-by-channels.

[~,~,ext] = fileparts(fileName);
ext = lower(ext);

switch ext
    case '.mat'
        s = load(fileName);
        if isfield(s,'amplifier_data')
            data = double(s.amplifier_data);
        elseif isfield(s,'raw_data')
            data = double(s.raw_data);
        else
            error('MAT file lacks amplifier_data or raw_data.');
        end
        if isfield(s,'fs')
            fs = s.fs;
        elseif isfield(s,'Fs_final')
            fs = s.Fs_final;
        elseif isfield(s,'frequency_parameters') && isfield(s.frequency_parameters,'amplifier_sample_rate')
            fs = s.frequency_parameters.amplifier_sample_rate;
        else
            error('Sampling rate not found in MAT file.');
        end
    case '.rhd'
        [data, fs] = read_rhd_simple(fileName);
    otherwise
        error('Unsupported file extension: %s', ext);
end

if size(data,1) ~= 6 && size(data,2) == 6
    data = data.';
end
if size(data,1) ~= 6
    error('Data must contain exactly 6 channels.');
end
end

function [data, fs] = read_rhd_simple(fileName)
%READ_RHD_SIMPLE Minimal Intan RHD reader for amplifier data.
%   [DATA, FS] = READ_RHD_SIMPLE(FILENAME) returns DATA as channels x
%   samples and FS in Hz.

fid = fopen(fileName,'r','ieee-le');
if fid < 0
    error('Could not open %s',fileName);
end
cleanup = onCleanup(@()fclose(fid));
magic = fread(fid,1,'uint32');
if magic ~= hex2dec('c6912702')
    error('Not an Intan RHD file.');
end
verMain = fread(fid,1,'int16');
verSec  = fread(fid,1,'int16');
if verMain == 1
    sBlock = 60;
else
    sBlock = 128;
end
fs = fread(fid,1,'single');
% Skip 9 values after sample rate
fseek(fid, 4*9 + 2, 'cof');            % skip dsp_enabled and freq fields
fseek(fid, 2, 'cof');                   % notch filter mode
fseek(fid, 4*2, 'cof');                 % impedance test freq
for k=1:3
    len = fread(fid,1,'uint32');
    fseek(fid,len,'cof');
end
if (verMain == 1 && verSec >= 1) || (verMain > 1)
    fseek(fid,2,'cof');                 % temp sensor channels
end
if (verMain == 1 && verSec >= 3) || (verMain > 1)
    fseek(fid,2,'cof');                 % board mode
end
if verMain > 1
    len = fread(fid,1,'uint32');
    fseek(fid,len,'cof');               % reference channel name
end
% Skip remainder of frequency_parameters structure (already read)
% Build channel structures
num_groups = fread(fid,1,'int16');
amp_ch = 0;
for g = 1:num_groups
    len = fread(fid,1,'uint32'); fseek(fid,len,'cof'); % group name
    len = fread(fid,1,'uint32'); fseek(fid,len,'cof'); % prefix
    fseek(fid,2,'cof');                              % enabled
    num_ch = fread(fid,1,'int16');
    num_amp = fread(fid,1,'int16'); %#ok<NASGU>
    for c = 1:num_ch
        len = fread(fid,1,'uint32'); fseek(fid,len,'cof'); % native name
        len = fread(fid,1,'uint32'); fseek(fid,len,'cof'); % custom name
        fseek(fid,2*7,'cof');                             % various int16 fields
        fseek(fid,4*2,'cof');                             % impedance
        enabled = fread(fid,1,'int16');
        if enabled
            type = fread(fid,1,'int16');
            if type==0
                amp_ch = amp_ch + 1;
            end
        else
            fseek(fid,2,'cof'); % skip unused type field
        end
    end
end
num_amp_channels = amp_ch;
bytes_per_block = sBlock*4;                               % time stamp
bytes_per_block = bytes_per_block + sBlock*2*num_amp_channels;
bytes_remaining = ftell(fid);
fseek(fid,0,'eof');
filesize = ftell(fid);
bytes_remaining = filesize - bytes_remaining;
num_blocks = floor(bytes_remaining / bytes_per_block);
num_samples = num_blocks*sBlock;
fseek(fid, ftell(fid)-bytes_remaining,'bof');
if num_amp_channels ~= 6
    error('File must contain 6 amplifier channels.');
end
data = zeros(num_amp_channels, num_samples, 'single');
for b = 1:num_blocks
    fseek(fid, sBlock*4, 'cof');
    raw = fread(fid,[sBlock,num_amp_channels],'uint16');
    data(:,(b-1)*sBlock+(1:sBlock)) = raw.';
end
data = 0.195*(double(data) - 32768);
end

function [eeg, t, fsOut] = preprocess(data, fsIn)
%PREPROCESS Filter and downsample EEG.
%   [EEG, T, FSOUT] = PREPROCESS(DATA, FSIN) returns re-referenced and
%   filtered EEG sampled at 1 kHz.

parietal = data(5,:);
active   = data(1:5,:) - parietal;
active(5,:) = 0;

hp = designfilt('highpassiir','FilterOrder',4, 'HalfPowerFrequency',1,'SampleRate',fsIn);
lp = designfilt('lowpassiir','FilterOrder',4, 'HalfPowerFrequency',250,'SampleRate',fsIn);
notch = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',58, 'HalfPowerFrequency2',62, 'SampleRate',fsIn);

eq = filtfilt(hp, active.');
eq = filtfilt(lp, eq); eq = filtfilt(notch, eq); eq = eq.';

fsOut = 1000;
if fsIn ~= fsOut
    eq = resample(eq.', fsOut, fsIn).';
end

t = (0:size(eq,2)-1).'/fsOut;

eeg = eq.'; %#ok<NASGU>
end

function bursts = detectBursts(eeg, fs)
%DETECTBURSTS Identify tic bursts in band 4-10 Hz.

bp = designfilt('bandpassiir','FilterOrder',4, 'HalfPowerFrequency1',4, ...
    'HalfPowerFrequency2',10,'SampleRate',fs);
filtered = filtfilt(bp, eeg);
env = abs(hilbert(filtered));
thr = mean(env) + 3*std(env);
mask = any(env(:,1:4) > thr(1:4),2);
minS = round(0.2*fs);
d = diff([0; mask; 0]);
startIdx = find(d==1); endIdx = find(d==-1)-1;
valid = (endIdx-startIdx+1)>=minS;
startIdx = startIdx(valid); endIdx = endIdx(valid);
N = numel(startIdx);
labels = ["R\_M2" "L\_M2" "R\_M1" "L\_M1"];
PeakTime = zeros(N,1); PeakPower = zeros(N,1); Dom = strings(N,1);
for k=1:N
    seg = startIdx(k):endIdx(k);
    [pk,idx] = max(env(seg,1:4),[],'all','linear');
    PeakPower(k) = pk.^2; %#ok<*NASGU>
    [idxT,chan] = ind2sub(size(env(:,1:4)),idx+seg(1)-1);
    PeakTime(k) = idxT/fs;
    Dom(k) = labels(chan);
end
bursts = table(startIdx/fs, endIdx/fs, PeakTime, PeakPower, Dom, ...
    'VariableNames',{'Onset','Offset','PeakTime','PeakPower','DominantChannel'});
end

function [coherence, phaseLag] = computeCoherence(eeg, bursts, fs)
%COMPUTECOHERENCE Coherence and phase metrics.

nCh = size(eeg,2);
coherence.burst = zeros(nCh,nCh);
coherence.baseline = zeros(nCh,nCh);

segments = round([bursts.Onset bursts.Offset]*fs);
mask = false(size(eeg,1),1);
for k=1:size(segments,1)
    mask(segments(k,1):segments(k,2)) = true;
end

pairs = nchoosek(1:nCh,2);
for i=1:size(pairs,1)
    x = eeg(:,pairs(i,1));
    y = eeg(:,pairs(i,2));
    [c,f] = mscohere(x,y,hamming(fs),fs/2,fs,fs);
    idx = f>=4 & f<=10;
    coherence.baseline(pairs(i,1),pairs(i,2)) = mean(c(idx));
    coherence.baseline(pairs(i,2),pairs(i,1)) = coherence.baseline(pairs(i,1),pairs(i,2));
    if any(mask)
        [c2,~] = mscohere(x(mask),y(mask),hamming(fs),fs/2,fs,fs);
        coherence.burst(pairs(i,1),pairs(i,2)) = mean(c2(idx));
        coherence.burst(pairs(i,2),pairs(i,1)) = coherence.burst(pairs(i,1),pairs(i,2));
    end
end

phaseLag.M1 = phase_at_peaks(eeg(:,3), eeg(:,4), bursts.PeakTime, fs);
phaseLag.M2 = phase_at_peaks(eeg(:,1), eeg(:,2), bursts.PeakTime, fs);
end

function lag = phase_at_peaks(chA, chB, times, fs)
%PHASE_AT_PEAKS Instantaneous phase difference.

bp = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',4, ...
    'HalfPowerFrequency2',10,'SampleRate',fs);
ha = hilbert(filtfilt(bp,chA));
hb = hilbert(filtfilt(bp,chB));
lag = zeros(numel(times),1);
for k=1:numel(times)
    idx = round(times(k)*fs);
    if idx>0 && idx<=length(ha)
        lag(k) = angle(ha(idx)/hb(idx));
    end
end
end

function bandPower = computeBandPower(eeg, fs)
%COMPUTEBANDPOWER Welch spectral power in classical bands.

bands = [4 8; 15 30; 30 80];
labels = {'Theta','Beta','Gamma'};
P = zeros(numel(labels), size(eeg,2));
for ch = 1:size(eeg,2)
    [pxx,f] = pwelch(eeg(:,ch), hamming(2*fs), [], [], fs);
    for b = 1:numel(labels)
        idx = f>=bands(b,1) & f<=bands(b,2);
        P(b,ch) = trapz(f(idx), pxx(idx));
    end
end
bandPower = struct('Theta',P(1,:), 'Beta',P(2,:), 'Gamma',P(3,:));
end

function plotSummary(t, eeg, bursts, coherence, fs)
%PLOTSUMMARY Generate figure summarizing analysis.

figure('Color','w','Position',[100 100 1200 800]);
chLabels = {'R M2','L M2','R M1','L M1','Parietal'};

subplot(3,2,1);
plot(t,eeg); xlabel('Time (s)'); ylabel('\muV'); title('Filtered EEG'); legend(chLabels);

subplot(3,2,3);
if ~isempty(bursts)
    win = round(fs*[-1 2]);
    [f,~,logS] = computeLogSpectrogram(eeg(:,1), fs, 256, 200, 80);
    img = zeros(length(f), diff(win)+1);
    for b = 1:height(bursts)
        idx = round(bursts.PeakTime(b)*fs) + win(1):round(bursts.PeakTime(b)*fs) + win(2);
        idx(idx<1 | idx>size(eeg,1)) = [];
        [~,~,tmp] = computeLogSpectrogram(eeg(idx,1), fs, 256, 200, 80);
        img(:,1:size(tmp,2)) = img(:,1:size(tmp,2)) + tmp;
    end
    img = img / height(bursts);
    imagesc((-1:2),f,img); axis xy; xlabel('Time (s)'); ylabel('Hz'); title('Burst Spectrogram');
else
    title('Burst Spectrogram - none'); axis off
end

subplot(3,2,5);
imagesc(coherence.burst); colorbar; caxis([0 1]); axis equal tight;
set(gca,'XTick',1:5,'XTickLabel',chLabels,'YTick',1:5,'YTickLabel',chLabels);
title('Burst Coherence');

subplot(3,2,6);
imagesc(coherence.baseline); colorbar; caxis([0 1]); axis equal tight;
set(gca,'XTick',1:5,'XTickLabel',chLabels,'YTick',1:5,'YTickLabel',chLabels);
title('Baseline Coherence');

dateStr = datestr(now,'yyyymmdd');
saveas(gcf, ['summary_' dateStr '.png']);
end

function [fOut,tOut,logSpec] = computeLogSpectrogram(sig, Fs_in, NFFT_val, NOvrlp, fCut)
%COMPUTELOGSPECTROGRAM Legacy log spectrogram calculation.

window_vec = hamming(NFFT_val);
[~,fAll,tOut,P_psd] = spectrogram(double(sig), window_vec, NOvrlp, NFFT_val, Fs_in, 'psd');
keep = fAll <= fCut;
fOut = fAll(keep);
P_subset = P_psd(keep,:);
min_psd_val = 1e-12;
P_subset(~isfinite(P_subset) | P_subset <= 0) = min_psd_val;
logSpec = log10(sqrt(P_subset));
end
%% --- END OF FILE ---
