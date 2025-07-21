% ESPRIT: eig() vs svd() DOA estimation comparison
clear; clc;

% Parameters
N = 8;                   % Number of sensors
% K = 2;                   % Number of sources fixed
T = 2;                 % Snapshots
true_angles = [-20 30];  % True DOAs in degrees
d = 0.5;                 % Element spacing (in wavelengths)
lambda = 1;              % Wavelength
k = 2*pi/lambda;         % Wavenumber
SNR_dB_list = -10:5:20;  % Range of SNRs to test

% Steering vector function
steering_vector = @(theta) exp(1j * k * d * (0:N-1).' * sind(theta));

% Initialize error arrays
err_eig = NaN(size(SNR_dB_list));
err_svd = NaN(size(SNR_dB_list));


% Main loop over SNR values
for i = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(i);

    % --- Signal generation ---
    A = zeros(N, K);
    for j = 1:K
        A(:, j) = steering_vector(true_angles(j));
    end
    % S = randn(K, T) + 1j * randn(K, T);      % Source signals
    % uncorrelated

    % Simulate correlated sources
    base_signal = randn(1, T) + 1j * randn(1, T);
    S = zeros(K, T);
    S(1, :) = base_signal;
    S(2, :) = base_signal + 0.05 * (randn(1, T) + 1j * randn(1, T));  % Highly correlated with small perturbation
    % You can add more sources similarly if K > 2


    noise = randn(N, T) + 1j * randn(N, T);  % Noise
    signal_power = mean(abs(S(:)).^2);
    noise_power = signal_power / (10^(SNR_dB/10));
    X = A*S + sqrt(noise_power)*noise;       % Received data

    %% --- EIG-based ESPRIT ---
    Rxx = (X * X') / T;
    [V, D] = eig(Rxx);
    [~, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    Vs = V(:, 1:K);  % Signal subspace

    % Solve shift-invariance equation
    Phi_eig = (Vs(1:end-1, :) \ Vs(2:end, :));
    eigvals_eig = eig(Phi_eig);
    doas_eig = real(asind(angle(eigvals_eig) / (2*pi*d)));
    doas_eig = doas_eig(isfinite(doas_eig));

    % Check and assign error if valid
    if numel(doas_eig) >= K
        doas_eig = sort(doas_eig(1:K));
        err_eig(i) = mean(abs(sort(true_angles(:)) - doas_eig(:)));
    else
        err_eig(i) = NaN;
    end

    %% --- SVD-based ESPRIT ---
    [U, ~, ~] = svd(X, 'econ');
    Us = U(:, 1:K);

    Phi_svd = (Us(1:end-1, :) \ Us(2:end, :));
    eigvals_svd = eig(Phi_svd);
    doas_svd = real(asind(angle(eigvals_svd) / (2*pi*d)));
    doas_svd = doas_svd(isfinite(doas_svd));

    if numel(doas_svd) >= K
        doas_svd = sort(doas_svd(1:K));
        err_svd(i) = mean(abs(sort(true_angles(:)) - doas_svd(:)));
    else
        err_svd(i) = NaN;
    end
end

%% --- Plotting results ---
validIdx = ~isnan(err_eig) & ~isnan(err_svd);

figure;
plot(SNR_dB_list(validIdx), err_eig(validIdx), 'r-o', 'LineWidth', 1.5); hold on;
plot(SNR_dB_list(validIdx), err_svd(validIdx), 'b--s', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('Mean DOA Estimation Error (degrees)');
title('ESPRIT DOA Estimation: eig() vs svd()');
legend('EIG-based ESPRIT', 'SVD-based ESPRIT', 'Location', 'northeast');
grid on;
