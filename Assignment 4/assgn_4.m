results = cell(1, 10);  % Store output from each run
times = cell(1, 10);    % Store time vector

for i = 1:10
    disp(['Running iteration ', num2str(i)]);

    % Run model
    sim('flann_architecture', 'StopTime', '0.0001668');

    % Extract timeseries data
    y = squeeze(ans.simout.Data);
    t = ans.simout.Time;

    % Save for this iteration
    results{i} = y;
    times{i} = t;

    % Display quick look
    disp(['First 5 output values of run ', num2str(i), ':']);
    disp(y(1:5));
end
%% Convert cell array to matrix
% Convert the 1x10 cell (each cell: 6005x1) to a 10x6005 matrix
% Preallocate
err_mat = zeros(10, 6005);

% Fill each row from the cell
for i = 1:10
    err_mat(i, :) = results{i}';  % Transpose column vector to row
end

%% Smoothing
disp('Please Wait! Smoothing Operation is Going On...');
length_of_smoothing_filter = 200;
smoothing_filter_coeff = (1 / length_of_smoothing_filter) * ones(1, length_of_smoothing_filter);

err_smooth = zeros(size(err_mat));
for i = 1:10
    err_smooth(i,:) = filter(smoothing_filter_coeff, 1, err_mat(i,:));
end
%% Plot Learning Curve
figure;
plot(10 * log10(mean(err_smooth))); 
xlabel('Iterations');
ylabel('MSE (dB)');
title('Learning Curve');
grid on;

%% Average MSE over last 1000 iterations
fln_mse = 10 * log10(mean(mean(err_mat(:, end-999:end))));
fprintf('Average MSE Value over the last 1000 iterations is %f dB\n', fln_mse);