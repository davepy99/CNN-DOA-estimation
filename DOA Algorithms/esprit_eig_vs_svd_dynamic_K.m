clear; clc;

% Parameters
N = 8;                   % Number of sensors
T = 5;                  % Snapshots (small T to stress the system)
true_angles = [-20 30];  % True DOAs (degrees)
d = 0.5;                 % Element spacing (wavelength units)
lambda = 1;              % Wavelength
k = 2*pi/lambda;
SNR_dB_list = -10:5:20;

% Steering vector function
steering_vector = @(theta) exp(1j * k * d * (0:N-1).' * sind(theta));

% Initialize errors
err_eig = NaN(size(SNR_dB_list));
err_svd = NaN(size(SNR_dB_list));
estimated_K = zeros(size(SNR_dB_list));

% Main loop
for i = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(i);

    % === Generate correlated sources ===
    base_signal = randn(1, T) + 1j * randn(1, T);
    S = zeros(2, T);
    S(1, :) = base_signal;
    S(2, :) = base_signal + 0.05 * (randn(1, T) + 1j * randn(1, T));  % Highly correlated

    A = zeros(N, 2);
    for j = 1:2
        A(:, j) = steering_vector(true_angles(j));
    end

    % Add noise
    noise = randn(N, T) + 1j * randn(N, T);
    signal_power = mean(abs(S(:)).^2);
    noise_power = signal_power / (10^(SNR_dB/10));
    X = A*S + sqrt(noise_power)*noise;

    % === SVD: Estimate K dynamically ===
    [U, S_svd, ~] = svd(X, 'econ');
    sing_vals = diag(S_svd);
    threshold = 0.1 * max(sing_vals);  % Heuristic threshold
    K_est = sum(sing_vals > threshold);
    K_est = min(K_est, N-1);  % Safety check
    estimated_K(i) = K_est;

    if K_est < 1
        continue;  % Skip invalid estimation
    end

    %% === SVD-based ESPRIT ===
    Us = U(:, 1:K_est);
    Phi_svd = (Us(1:end-1, :) \ Us(2:end, :));
    eigvals_svd = eig(Phi_svd);
    doas_svd = real(asind(angle(eigvals_svd) / (2*pi*d)));
    doas_svd = doas_svd(isfinite(doas_svd));

    if numel(doas_svd) >= 2
        doas_svd = sort(doas_svd(1:2));
        err_svd(i) = mean(abs(sort(true_angles(:)) - doas_svd(:)));
    end

    %% === EIG-based ESPRIT ===
    Rxx = (X * X') / T;
    [V, D] = eig(Rxx);
    [~, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    Vs = V(:, 1:K_est);

    Phi_eig = (Vs(1:end-1, :) \ Vs(2:end, :));
    eigvals_eig = eig(Phi_eig);
    doas_eig = real(asind(angle(eigvals_eig) / (2*pi*d)));
    doas_eig = doas_eig(isfinite(doas_eig));

    if numel(doas_eig) >= 2
        doas_eig = sort(doas_eig(1:2));
        err_eig(i) = mean(abs(sort(true_angles(:)) - doas_eig(:)));
    end
end

%% === Plot DOA estimation error ===
figure;
validIdx = ~isnan(err_eig) & ~isnan(err_svd);
plot(SNR_dB_list(validIdx), err_eig(validIdx), 'r-o', 'LineWidth', 1.5); hold on;
plot(SNR_dB_list(validIdx), err_svd(validIdx), 'b--s', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('Mean DOA Estimation Error (degrees)');
title('ESPRIT with Dynamic K: EIG vs SVD');
legend('EIG-based ESPRIT', 'SVD-based ESPRIT', 'Location', 'northeast');
grid on;

%% === Plot estimated K ===
figure;
plot(SNR_dB_list, estimated_K, 'k-x', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('Estimated Number of Sources');
title('Dynamic Rank Estimation from SVD');
grid on;

%% === Optional: Plot singular value spectrum at a specific SNR ===
% You can run this once more for a specific SNR index:
snr_index = find(SNR_dB_list == 0, 1);
if ~isempty(snr_index)
    % Re-run SVD for SNR = 0
    SNR_dB = SNR_dB_list(snr_index);
    base_signal = randn(1, T) + 1j * randn(1, T);
    S = [base_signal; base_signal + 0.05*(randn(1,T)+1j*randn(1,T))];
    A = [steering_vector(true_angles(1)), steering_vector(true_angles(2))];
    noise = randn(N, T) + 1j * randn(N, T);
    signal_power = mean(abs(S(:)).^2);
    noise_power = signal_power / (10^(SNR_dB/10));
    X = A*S + sqrt(noise_power)*noise;

    [~, S_svd_demo, ~] = svd(X, 'econ');
    figure;
    semilogy(diag(S_svd_demo), 'ko-', 'LineWidth', 1.5);
    xlabel('Index');
    ylabel('Singular Value (log scale)');
    title('Singular Value Spectrum at SNR = 0 dB');
    grid on;
end
