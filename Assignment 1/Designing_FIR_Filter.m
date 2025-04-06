%% 1. SIGNAL LOADING AND PREPROCESSING
% Load both original and recorded signals
[original, fs_orig] = audioread('MLS_excitation_FINAL.wav'); % MLS Original File
[recorded, fs_rec] = audioread('/MATLAB Drive/MLSprat.wav');   % MLS recorded file from the phone
fs = fs_rec; % Use the recording sample rate for all calculations

% Display signal information
disp(['Original signal length: ', num2str(length(original)), ' samples (', num2str(length(original)/fs), ' seconds)']);
disp(['Recorded signal length: ', num2str(length(recorded)), ' samples (', num2str(length(recorded)/fs), ' seconds)']);

%% 2. PARAMETERS SETUP
% Define key parameters. These parameters are used to align the original MLS signal with the recorded MLS signal.
n = 17;                  
L = 2^n - 1;             % Length of the MLS sequence
numRuns = 4;             % Number of MLS runs
pre_delay = 3;           % Delay before the MLS signal starts (in seconds)
click_dur = 0.2;         % Duration of the initial click (in seconds)
post_click_space = 0.5;  % Pause after the initial click (in seconds)

%% 3. ORIGINAL MLS EXTRACTION
% Calculate the start of MLS in the original signal
mls_start_orig = round((pre_delay + click_dur + post_click_space) * fs) + 1;
mls_end_orig = min(mls_start_orig + L - 1, length(original));
mlsOriginal = original(mls_start_orig:mls_end_orig);

% Ensure that the number of samples equals L. This error may be due to a change in 
% the sampling rate during the creation or while reading of the original file.
if length(mlsOriginal) ~= L
    warning('MLS original length mismatch. Truncating or padding to the exact length.');
    if length(mlsOriginal) > L
        mlsOriginal = mlsOriginal(1:L);
    else
        mlsOriginal = [mlsOriginal; zeros(L - length(mlsOriginal), 1)];
    end
end

%% 4. CLICK DETECTION USING ENERGY DETECTION
disp('Searching for click using energy detection...');
% Normalize signals
original_norm = original / max(abs(original));
recorded_norm = recorded / max(abs(recorded));

% Parameters for click detection
click_start = round(pre_delay * fs) + 1;
click_end = click_start + round(click_dur * fs) - 1;

% Calculate energy function (using short-time energy)
window_size = round(0.01 * fs); % 10ms window
energy_recorded = zeros(length(recorded_norm) - window_size, 1);
for i = 1:length(energy_recorded)
    energy_recorded(i) = sum(recorded_norm(i:i+window_size-1).^2);
end

% Find peaks in energy function
threshold = 0.008 * max(energy_recorded); % Adaptive threshold according to the distance of the two phones 
[peaks, peak_locations] = findpeaks(energy_recorded, 'MinPeakHeight', threshold);

% If no peaks found, try with lower threshold
if isempty(peaks)
    threshold = 0.002 * max(energy_recorded);
    [peaks, peak_locations] = findpeaks(energy_recorded, 'MinPeakHeight', threshold);
end

% Sort peaks by position (time)
[sorted_locations, sorted_idx] = sort(peak_locations, 'ascend');
sorted_peaks = peaks(sorted_idx);

% Display top candidates for manual inspection
disp('Top 5 peak candidates:');
num_candidates = min(5, length(sorted_locations));
for i = 1:num_candidates
    disp(['Candidate ', num2str(i), ': Position=', num2str(sorted_locations(i))]);
end

% Use the first peak (based on temporal order)
if ~isempty(sorted_locations)
    click_position = sorted_locations(1);
    disp(['Click detected at sample: ', num2str(click_position)]);
    
    % Calculate MLS start position using the known delay after the click
    % In original signal: MLS starts post_click_space seconds after click end
    mls_start_rec = click_position + round((click_dur + post_click_space) * fs);
    disp(['MLS start calculated at sample: ', num2str(mls_start_rec)]);
else
    error('No click detected. Try adjusting the threshold parameters.');
end

%% 5. MLS EXTRACTION
% Ensure MLS start is within valid range
if mls_start_rec <= 0
    mls_start_rec = 1;
    warning('Calculated MLS start position is before recording start. Using beginning of recording.');
elseif mls_start_rec > length(recorded)
    error('Calculated MLS start position exceeds recording length.');
end

% Calculate the end position of the MLS extraction
required_length = L * numRuns; % Length for all MLS repetitions
mls_end_rec = min(mls_start_rec + required_length - 1, length(recorded));

% Extract the MLS portion from the recorded signal (removing click and pre-delay)
mls_recorded = recorded(mls_start_rec:mls_end_rec);
disp(['Extracted MLS length: ', num2str(length(mls_recorded)), ' samples']);
disp(['Number of complete MLS runs extracted: ', num2str(floor(length(mls_recorded)/L))]);

% Check if we got enough samples for all repetitions
if length(mls_recorded) < L
    error('Not enough MLS data extracted for even one repetition');
elseif length(mls_recorded) < L*numRuns
    warning(['Only extracted ', num2str(floor(length(mls_recorded)/L)), ' complete MLS repetitions']);
    actual_runs = floor(length(mls_recorded)/L);
else
    actual_runs = numRuns;
end

%% 6. VISUALIZATION - SIMPLIFIED TO JUST TWO SUBPLOTS
% Plot original and recorded signals with MLS start markers
figure(1);
subplot(2,1,1);
plot(original_norm);
hold on;
plot([click_start, click_start], [-1, 1], 'r', 'LineWidth', 2);
plot([mls_start_orig, mls_start_orig], [-1, 1], 'g', 'LineWidth', 2);
title('Original Signal with Click and MLS Start');
xlabel('Sample');
legend('Signal', 'Click Start', 'MLS Start');
grid on;

subplot(2,1,2);
plot(recorded_norm);
hold on;
plot([click_position, click_position], [-1, 1], 'r', 'LineWidth', 2);
plot([mls_start_rec, mls_start_rec], [-1, 1], 'g', 'LineWidth', 2);
title('Recorded Signal with Click and MLS Start');
xlabel('Sample');
legend('Signal', 'Click Position', 'MLS Start');
grid on;


%% HELPER FUNCTION
function ir = fft_deconvolve(recorded, original)
    % Function for MLS deconvolution using FFT
    N = length(original);
    X = fft(recorded);
    Y = fft(original);
    H = X ./ Y;
    ir = real(ifft(H));

    % Circular shift to get causal IR
    [~, max_idx] = max(abs(ir));
    ir = circshift(ir, -max_idx+1);
end
%% 7. IMPULSE RESPONSE CALCULATION
% MLS deconvolution using circular cross-correlation
ir = zeros(L, 1);

% Process each MLS repetition and average to improve SNR
for i = 1:numRuns
    start_idx = (i-1)*L + 1;
    end_idx = min(i*L, length(mls_recorded));
    
    if end_idx - start_idx + 1 < L
        % Pad with zeros if this segment is too short
        segment = [mls_recorded(start_idx:end_idx); zeros(L - (end_idx - start_idx + 1), 1)];
    else
        segment = mls_recorded(start_idx:start_idx+L-1);
    end
    
    segment_ir = fft_deconvolve(segment, mlsOriginal);
    ir = ir + segment_ir;
end

% Average the impulse responses
ir = ir / numRuns;

%% 8. IMPULSE RESPONSE VISUALIZATION AND ANALYSIS
% Plot impulse response
figure(3);
t_ir = (0:length(ir)-1)/fs;
plot(t_ir, ir);
title('Room Impulse Response');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% 9. INVERSE FILTER CREATION
% Create the frequency-domain inverse filter
fft_size = 2048; % Adjust based on IR length
H = fft(ir, fft_size);

% Plot frequency response (optional)
figure(4);
f = (0:fft_size/2)*fs/fft_size;
magnitude = 20*log10(abs(H(1:fft_size/2+1)));
plot(f, magnitude);
title('Frequency Response of Room');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

% Regularization to prevent division by small values
lambda = 0.01; % Adjust based on your SNR
H_inverse = conj(H) ./ (abs(H).^2 + lambda);

% Create the frequency domain filter object
fdf = dsp.FrequencyDomainFIRFilter('NumeratorDomain', 'Frequency', ...
                                  'Numerator', H_inverse, ...
                                  'PartitionForReducedLatency', false);
%% 10. FILTER ANALYSIS
% Create the frequency-domain inverse filter (from your code)
fft_size = 2048; % Adjust based on IR length
H = fft(ir, fft_size);

% Regularization to prevent division by small values
lambda = 0.01; % Adjust based on your SNR
H_inverse = conj(H) ./ (abs(H).^2 + lambda);

% Create the frequency domain filter object
fdf = dsp.FrequencyDomainFIRFilter('NumeratorDomain', 'Frequency', ...
                                  'Numerator', H_inverse, ...
                                  'PartitionForReducedLatency', false);

%% 11. FILTER ANALYSIS AND VISUALIZATION
freq_coeffs = fdf.Numerator;% Here we are extracting frequency filter coefficients from the Filter 
time_coeffs = ifft(freq_coeffs);% Convert frequency domain coefficients to time domain
time_coeffs = real(time_coeffs); 
firFilt = dsp.FIRFilter('Numerator', time_coeffs(:).');
figure(5);
zplane(firFilt.Numerator, 1);% Create pole-zero plot (all poles at origin for FIR filter)
title('Pole-Zero Plot of Room Inverse Filter');
grid on;

figure(6);
subplot(3,1,1);
[h, w] = freqz(fdf.Numerator, 1, 1024, fs);% Plot Magnitude Response
plot(w, 20*log10(abs(h)));
title('Magnitude Response of Inverse Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

subplot(3,1,2);
phase = unwrap(angle(h))*180/pi;% Plot Phase Response
plot(w, phase);
title('Phase Response of Inverse Filter');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
grid on;

subplot(3,1,3);
gd = grpdelay(time_coeffs, 1, 1024, fs);% Plot Group Delay
plot(w, gd);
title('Group Delay of Inverse Filter');
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
grid on;

% Displaying filter properties and stability analysis
disp('---- Filter Analysis ----');
disp(['Filter Type: Frequency Domain FIR']);
disp(['Numerator Domain: ' fdf.NumeratorDomain]);
disp(['Filter Method: ' fdf.Method]);
disp(['FFT Size: ' num2str(fft_size)]);
disp(['Partition For Reduced Latency: ' num2str(fdf.PartitionForReducedLatency)]);

% Checking filter stability based on coefficient magnitudes
max_coeff = max(abs(time_coeffs));
disp(['Maximum coefficient magnitude: ' num2str(max_coeff)]);
if max_coeff > 10
    warning('Filter may have stability issues due to large coefficients');
else
    disp('Filter appears stable based on coefficient magnitudes');
end

% Calculating filter energy
energy = sum(abs(time_coeffs).^2);
disp(['Filter energy: ' num2str(energy)]);
%% 10. SPEECH COMPENSATION
% Load speech file
[speech, fs_speech] = audioread('/MATLAB Drive/Modi Ji Sanskrit Speaking Moment_Solo.wav');
% Apply inverse filter to compensate for room acoustics
disp('Applying inverse filter to speech...');
compensated_speech = fdf(speech);

% Normalize the output to avoid clipping
max_val = max(abs(compensated_speech));
if max_val > 0.99
    compensated_speech = compensated_speech / max_val * 0.95;
    disp('Normalized compensated speech to avoid clipping');
end

%% 11. SAVE AND COMPARE RESULTS
% Save compensated speech for playback
output_file = '/MATLAB Drive/Modi Ji Sanskrit Speaking Moment_Compensated_final_prats.wav';
audiowrite(output_file, compensated_speech, fs);
disp(['Compensated speech saved to: ', output_file]);


