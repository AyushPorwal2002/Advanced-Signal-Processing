%% 1. GENERATE MLS SIGNAL
clear all; close all; clc;

% Parameters
fs = 44100;             % Sampling rate (Hz)
n = 17;                 % Order of the MLS sequence
L = 2^n - 1;            % Length of MLS sequence
numRuns = 4;            % Number of repetitions
pre_delay = 3;          % 3-second lead silence
click_dur = 0.2;        % 200ms double click Detection Zone 
post_click_space = 0.5; % 500ms Post click zeroes
level = -5;             % Excitation level in dB

mlsOriginal = mls(L, 'ExcitationLevel', level); % Create MLS sequence of order n

% Build synchronization header so that there is no delay while recording and playing the signal (human error)
header = [zeros(pre_delay*fs,1);             % Leading zeroes
          ones(click_dur*fs,1)*0.8;          % 200ms Double Click Detection Zone
          zeros(post_click_space*fs,1)];     % Post-click zeroes

% Full excitation signal
excitation = [header; 
              repmat(mlsOriginal, numRuns, 1); 
              zeros(L + 1,1)];

% Save as WAV file
audiowrite('MLS_excitation_FINAL.wav', excitation, fs);
%% Plotting Section
figure('Color','white','Position',[100 100 900 600])

% Time domain plot
subplot(2,1,1)
t = (0:length(excitation)-1)/fs; % Time vector in seconds
plot(t, excitation)
title('MLS Excitation Signal - Time Domain')
xlabel('Time (seconds)')
ylabel('Amplitude')
xlim([0 t(end)])
grid on

% Highlight different sections
hold on
plot([pre_delay pre_delay],[-1 1],'r--') % Start of click
plot([pre_delay+click_dur pre_delay+click_dur],[-1 1],'r--') % End of click
hold off
legend('Signal','Section Markers')

% Frequency domain plot
subplot(2,1,2)
mls_part = repmat(mlsOriginal, numRuns, 1);
N_mls = length(mls_part);
f_mls = fs*(0:(N_mls/2))/N_mls;
Y_mls = fft(mls_part);
P_mls = abs(Y_mls/N_mls).^2;
PdB_mls = 10*log10(P_mls(1:N_mls/2+1));

% Plot the corrected spectrum
plot(f_mls, PdB_mls)
title('Frequency Spectrum')
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
set(gca,'XScale','log')
grid on
xlim([20 20000]) % Human hearing range
