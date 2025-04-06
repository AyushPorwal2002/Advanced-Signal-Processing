%% Load audio files
ref_file = '/MATLAB Drive/Modi Ji Sanskrit Speaking Moment_Solo.wav';
deg_file_before = '/MATLAB Drive/modipratsoriginsl.wav';
deg_file_after = '/MATLAB Drive/finalpratsfiltermodi.wav';

[ref_data, ref_fs] = audioread(ref_file);
[deg_data_before, fs_before] = audioread(deg_file_before);
[deg_data_after, fs_after] = audioread(deg_file_after);

%% Time axis generation
ref_time = (0:length(ref_data)-1) / ref_fs;
before_time = (0:length(deg_data_before)-1) / fs_before;
after_time = (0:length(deg_data_after)-1) / fs_after;

%% Plot Original Signal
figure;
subplot(2,1,1);
plot(ref_time, ref_data, 'b');
legend('Original Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Original Signal');
grid on;

%% Plot Recorded Signal Before Filtering
subplot(2,1,2);
plot(before_time, deg_data_before, 'r');
legend('Recorded Signal (Before Filtering)');
xlabel('Time (s)'); ylabel('Amplitude');
title('Recorded Signal Before Filtering');
grid on;

%% Plot Original Signal Again for After Filtering Comparison
figure;
subplot(2,1,1);
plot(ref_time, ref_data, 'b');
legend('Original Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Original Signal');
grid on;

%% Plot Filtered Signal
subplot(2,1,2);
plot(after_time, deg_data_after, 'g');
legend('Filtered Signal');
xlabel('Time (s)'); ylabel('Amplitude');
title('Filtered Signal');
grid on;

%% Frequency Spectrum of Original, Recorded, and Filtered Signal
figure;
subplot(3,1,1);
freq_ref = abs(fft(ref_data));
freq_axis_ref = linspace(0, ref_fs/2, length(freq_ref)/2);
plot(freq_axis_ref, freq_ref(1:length(freq_axis_ref)), 'b');
legend('Original Signal Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Spectrum of Original Signal');
grid on;

subplot(3,1,2);
freq_deg_before = abs(fft(deg_data_before));
freq_axis_before = linspace(0, fs_before/2, length(freq_deg_before)/2);
plot(freq_axis_before, freq_deg_before(1:length(freq_axis_before)), 'r');
legend('Recorded Signal Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Spectrum of Recorded Signal Before Filtering');
grid on;

subplot(3,1,3);
freq_deg_after = abs(fft(deg_data_after));
freq_axis_after = linspace(0, fs_after/2, length(freq_deg_after)/2);
plot(freq_axis_after, freq_deg_after(1:length(freq_axis_after)), 'g');
legend('Filtered Signal Spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Spectrum of Filtered Signal');
grid on;
