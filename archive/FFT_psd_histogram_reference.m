% Project : Tic treatment via ultrasound stimulation
% Created by : Sanghyun Cho
% Email : daniel151@kaist.ac.kr
% Purpose : EEG analysis in Tic mouse model (pilot test)

%% initialization
clear; close all; clc;
%frequency & epoch
fs = 1000;
epoch = 6;

% Load the data
[file, path] = uigetfile('*.mat', 'Select a MATLAB data file', 'MultiSelect', 'off');
data = load(fullfile(path, file));
raw_eeg_data = data.amplifier_data(1,:);

% 파일 이름에서 .mat 확장자 제거
[~, name, ~] = fileparts(file);

%% raw data plot

time = (0:length(raw_eeg_data)-1) / fs; %time vector
figure()
plot(time, raw_eeg_data);
xlabel('시간 (초)');
xlim([0 round(length(raw_eeg_data)/fs)]); 
ylabel('EEG 신호');
ylim([-1500 1500]);
grid on;


%% Create spectrogram
figure()
spectrogram(raw_eeg_data, 300, 80, 300, fs, 'yaxis');
ylim([0 50]); 
 
%% freqeuncy parameter set up
cutoff_low_delta = 0.5;    % unit : Hz
cutoff_high_delta = 4;     % unit : Hz
cutoff_low_theta = 4;      % unit : Hz
cutoff_high_theta = 8;     % unit : Hz 
cutoff_low_alpha = 8;      % unit : Hz
cutoff_high_alpha = 12;    % unit : Hz
cutoff_low_beta = 12;      % unit : Hz
cutoff_high_beta = 30;     % unit : Hz
cutoff_low_gamma = 30;     % unit : Hz
cutoff_high_gamma = 100;   % unit : Hz

%% FFT analysis (devide by epoch & plot by frequency)
% delta, theta, alpha, beta, gamma가 시간에 따라 어떻게 변하는지 

data_total_duration = length(raw_eeg_data(1,:))/fs;
num_epoch = round(data_total_duration/epoch);
clear epoch_eeg; 
for i =1:num_epoch
    epoch_eeg(i,:) = raw_eeg_data(1+fs*epoch*(i-1):fs*epoch*i);
end
L = epoch*fs;
t = (0:epoch*fs-1)/fs;
n = 2^nextpow2(L);

ii=1;

for i = 1:num_epoch
    fft_eeg = fft(epoch_eeg(i,:),n);
    f_range = fs*(0:(n/2))/n;
    fft_eeg_mag = abs(fft_eeg/n);

    power_delta = 0;
    power_theta = 0;
    power_alpha = 0;
    power_beta = 0;
    power_gamma = 0;

    for j = 1:round(length(f_range)*101/(fs/2)) %FFT frequency component sweep for 0~100Hz
        fx = f_range(j);

        %delta frequency range value
        if fx>cutoff_low_delta && fx<cutoff_high_delta
            power_delta = power_delta + fft_eeg_mag(j);            
        end

        %theta frequency range value
        if fx>cutoff_low_theta && fx<cutoff_high_theta
            power_theta = power_theta + fft_eeg_mag(j);
        end

        %alpha frequency range value
        if fx>cutoff_low_alpha && fx<cutoff_high_alpha
            power_alpha = power_alpha + fft_eeg_mag(j);
        end

        %beta frequency range value
        if fx>cutoff_low_beta && fx<cutoff_high_beta
            power_beta = power_beta + fft_eeg_mag(j);
        end

         %gamma frequency range value
        if fx>cutoff_low_gamma && fx<cutoff_high_gamma
            power_gamma = power_gamma + fft_eeg_mag(j);
        end
    end

    power_epoch_value_delta(i) = power_delta;
    power_epoch_value_theta(i) = power_theta;
    power_epoch_value_alpha(i) = power_alpha;
    power_epoch_value_beta(i) = power_beta;
    power_epoch_value_gamma(i) = power_gamma;

end
freq_labels = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
power_values = [power_delta, power_theta, power_alpha, power_beta, power_gamma];

%% Power Spectral Density
figure()
bar(power_values);
set(gca, 'XTickLabel', freq_labels);
xlabel('Frequency Band');
ylabel('Power');
title('EEG Power Spectral Density');
grid on;

%% Plot data by frequency
figure()
plot(power_epoch_value_delta)
grid on;
figure()
plot(power_epoch_value_theta)
grid on;
figure()
plot(power_epoch_value_alpha)
grid on;
figure()
plot(power_epoch_value_beta)
grid on;
figure()
plot(power_epoch_value_gamma)
grid on;

%% FFT analysis (by frequency)
% FFT 수행
L_second = length(raw_eeg_data);           % 데이터 길이
n_second= 2^nextpow2(L_second);                  % 다음 2의 제곱수
fft_data_second = fft(raw_eeg_data, n_second);    % FFT 계산
f_second = fs * (0:(n_second/2)) / n_second;              % 주파수 벡터 생성

% FFT 결과의 크기 계산
fft_mag_second = abs(fft_data_second/n_second);          % 정규화된 FFT 크기
fft_mag_second = fft_mag_second(1:n_second/2+1);         % 0부터 Nyquist 주파수까지의 값만 사용
fft_mag_second(2:end-1) = 2 * fft_mag_second(2:end-1); % 양의 주파수 대역에 대해 크기 조정

% FFT 결과 시각화
figure;
plot(f_second, fft_mag_second);
xlim([0 100]); 
ylim([0 2]);
grid on;