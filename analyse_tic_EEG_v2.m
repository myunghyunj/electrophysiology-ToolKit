%% analyse_tic_EEG_v2.m
% Master script for epidural screw EEG analysis using six Intan inputs
% Requires MATLAB R2024a and Signal Processing Toolbox
%
% Channel mapping (Intan order):
% 1 - R_M2, 2 - L_M2, 3 - R_M1, 4 - L_M1, 5 - Parietal (reference), 6 - Ground
%
% The script loads .rhd or exported .mat files, re-references channels,
% performs filtering, burst detection and network analysis, then saves
% results and a summary figure.

%% Select recording file
[file,path] = uigetfile({'*.rhd;*.mat'},'Select Intan RHD or MAT file');
if isequal(file,0); error('No file selected'); end
fname = fullfile(path,file);

%% Load data
[raw,fs] = load_intan(fname);
if size(raw,1)~=6
    error('Expected 6 channels (5 active + ground). Found %d',size(raw,1));
end

%% Re-reference to parietal screw (channel 5)
rawRef = raw(1:5,:) - raw(5,:);
rawRef(5,:) = 0;

%% Pre-process
[eeg,t,fs_ds] = preprocess(rawRef,fs);

%% Detect bursts
bursts = detectBursts(eeg(:,1:4),fs_ds);

%% Network metrics
[bandPower,coherence,phaseLag] = computeCoherence(eeg,bursts,fs_ds);

%% Visualisation
plotSummary(t,rawRef,eeg,bursts,fs_ds,coherence);

%% Save results
save('results.mat','fs','t','eeg','bursts','bandPower','coherence','phaseLag');

%% ------- Local functions -------
function [data,fs] = load_intan(filename)
%LOAD_INTAN Load Intan RHD or MAT data.
%   [data,fs] = load_intan(filename)
ext = lower(extractAfter(filename,find(filename=='.',1,'last')));
if strcmp(ext,'rhd')
    S = read_Intan_RHD2000_file(filename);
    data = single(S.amplifier_data);
    fs = S.frequency_parameters.amplifier_sample_rate;
elseif strcmp(ext,'mat')
    S = load(filename);
    if isfield(S,'amplifier_data') && isfield(S,'fs')
        data = single(S.amplifier_data);
        fs = double(S.fs);
    else
        error('MAT file must contain variables amplifier_data and fs');
    end
else
    error('Unsupported file extension: %s',ext);
end
end

function [eeg,t,fsOut] = preprocess(x,fsIn)
%PREPROCESS Filter, notch and downsample signals.
fsOut = 1000;
[b,a] = butter(4,[1 250]/(fsIn/2),'bandpass');
notch = designfilt('bandstopiir','FilterOrder',2,...
    'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,...
    'DesignMethod','butter','SampleRate',fsIn);
xf = filtfilt(b,a,double(x).').';
xf = filtfilt(notch,double(xf).').';

eeg = resample(single(xf).',fsOut,fsIn).';
t = (0:size(eeg,2)-1)/fsOut;

eeg = eeg.'; % time x channels
end

function bursts = detectBursts(sig,fs)
%DETECTBURSTS Detect tic rhythm bursts using 4-10 Hz envelope.
[bpB,bpA] = butter(4,[4 10]/(fs/2),'bandpass');
filtSig = filtfilt(bpB,bpA,sig);
env = abs(hilbert(filtSig));
thr = mean(env) + 3*std(env);
mask = env>thr;
minLen = round(0.2*fs);

onset = []; offset = []; peak = []; chan = [];
for c=1:size(sig,2)
    m = mask(:,c);
    d = diff([0;m;0]);
    on = find(d==1); off = find(d==-1)-1;
    keep = (off-on+1)>=minLen;
    on = on(keep); off = off(keep);
    for k=1:numel(on)
        seg = env(on(k):off(k),:);
        [pk,idx] = max(seg(:));
        [~,dom] = ind2sub(size(seg),idx);
        onset(end+1) = (on(k)-1)/fs; %#ok<AGROW>
        offset(end+1) = (off(k)-1)/fs; %#ok<AGROW>
        peak(end+1) = pk; %#ok<AGROW>
        chan(end+1) = dom; %#ok<AGROW>
    end
end
bursts = table(onset',offset',peak',chan',...
    'VariableNames',{'onset','offset','peakPower','channel'});
end

function [bandPower,coherence,phaseLag] = computeCoherence(eeg,bursts,fs)
%COMPUTECOHERENCE Compute power, coherence and phase lag metrics.
ch = size(eeg,2);
% Spectral power
bands = [4 8;15 30;30 80];
P = zeros(3,ch);
for b=1:3
    [pxx,f] = pwelch(eeg,fs,[],[],fs);
    idx = f>=bands(b,1) & f<=bands(b,2);
    P(b,:) = trapz(f(idx),pxx(idx,:));
end
bandPower = P;

% Coherence during bursts and baseline
pairs = nchoosek(1:ch,2);
cohBurst = zeros(ch); cohBase = zeros(ch);
for p=1:size(pairs,1)
    i=pairs(p,1); j=pairs(p,2);
    cB = []; cR = [];
    for b=1:height(bursts)
        idx = round(bursts.onset(b)*fs):round(bursts.offset(b)*fs);
        base = idx - round(2*fs);
        base(base<1) = [];
        if numel(idx)<4, continue; end
        [c,f] = mscohere(eeg(idx,i),eeg(idx,j),[],[],[],fs);
        cB = [cB; c];
        [c,f] = mscohere(eeg(base,i),eeg(base,j),[],[],[],fs);
        cR = [cR; c];
    end
    fMask = f>=4 & f<=10;
    cohBurst(i,j) = mean(cB(:,fMask),'all');
    cohBase(i,j)  = mean(cR(:,fMask),'all');
    cohBurst(j,i)=cohBurst(i,j); cohBase(j,i)=cohBase(i,j);
end
coherence.burst = cohBurst; coherence.baseline = cohBase;

% Phase lag at burst peaks
bp = designfilt('bandpassiir','FilterOrder',4,
    'HalfPowerFrequency1',4,'HalfPowerFrequency2',10,
    'SampleRate',fs);
fSig = filtfilt(bp,eeg);
hSig = hilbert(fSig);
phaseLag.M1 = zeros(height(bursts),1);
phaseLag.M2 = zeros(height(bursts),1);
for b=1:height(bursts)
    idx = round((bursts.onset(b)+bursts.offset(b))*fs/2);
    if idx<1 || idx>size(eeg,1), continue; end
    phaseLag.M2(b) = angle(exp(1i*(angle(hSig(idx,1))-angle(hSig(idx,2)))));
    phaseLag.M1(b) = angle(exp(1i*(angle(hSig(idx,3))-angle(hSig(idx,4)))));
end
end

function plotSummary(t,rawRef,eeg,bursts,fs,coh)
%PLOTSUMMARY Generate multi-panel summary figure.
figure('Color','w','Position',[100 100 1200 800]);

subplot(3,2,1)
plot(t,rawRef'); axis tight
xlabel('Time (s)'); ylabel('\muV');
title('Re-referenced raw signals');

subplot(3,2,2)
plot(t,eeg); axis tight
xlabel('Time (s)'); ylabel('\muV');
title('Filtered & downsampled');

subplot(3,2,3)
win = round(3*fs);
for b=1:min(20,height(bursts))
    idx = round(bursts.onset(b)*fs)-fs:round(bursts.offset(b)*fs)+2*fs;
    idx(idx<1) = []; idx(idx>length(t))=[];
    [f,tt,ls] = computeLogSpectrogram(eeg(idx,1),fs,round(0.256*fs),round(0.9*0.256*fs),80);
    imagesc(tt-1,f,ls); hold on
end
axis xy; xlabel('Time (s)'); ylabel('Hz'); title('Burst aligned spectrogram');

subplot(3,2,5)
imagesc(coh.burst); axis square; caxis([0 1]); colorbar
xlabel('Channel'); ylabel('Channel'); title('Coherence during bursts');

subplot(3,2,6)
imagesc(coh.baseline); axis square; caxis([0 1]); colorbar
xlabel('Channel'); ylabel('Channel'); title('Baseline coherence');

sgtitle(datestr(now));
print(sprintf('summary_%s.png',datestr(now,'yyyymmdd')),'-dpng','-r150');
end

function [fOut,tOut,logSpec] = computeLogSpectrogram(sig,Fs,NFFT,NOv,fCut)
%COMPUTELOGSPECTROGRAM Helper for log spectrogram
w = hamming(NFFT);
[~,f,t,P] = spectrogram(double(sig),w,NOv,NFFT,Fs,'psd');
keep = f<=fCut; fOut=f(keep); P=P(keep,:);
P(P<=0)=1e-12; logSpec = log10(sqrt(P)); tOut=t;
end

%% --- END OF FILE ---
