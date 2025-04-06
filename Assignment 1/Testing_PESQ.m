
%% Calculate PESQ scores between reference and degraded audio files before applying filter 
ref_file = '/MATLAB Drive/Modi Ji Sanskrit Speaking Moment_Solo.wav';  
deg_file = '/MATLAB Drive/modipratsoriginsl.wav';    
 [ ref_data, ref_sampling_rate ] = audioread( ref_file ); 
 [ deg_data, deg_sampling_rate ] = audioread( deg_file );
scores = pesq(ref_file, deg_file);   % Calculation of PESQ by provided original and recorded file
fprintf('Wideband MOS-LQO before applying filter: %.3f\n', scores);
%% Calculate PESQ scores between reference and degraded audio files after applying filter
ref_file = '/MATLAB Drive/Modi Ji Sanskrit Speaking Moment_Solo.wav'; 
deg_file = '/MATLAB Drive/finalpratsfiltermodi.wav';  
 [ ref_data, ref_sampling_rate ] = audioread( ref_file ); 
 [ deg_data, deg_sampling_rate ] = audioread( deg_file );
scores = pesq(ref_file, deg_file); % % Calculation of PESQ by provided original and recorded file after filtering process
fprintf('Wideband MOS-LQO after applying filter: %.3f\n', scores);