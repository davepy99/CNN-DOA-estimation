clear; clc;

% Parameters
N = 8;                  % Number of sensors
T = 2;                 % Snapshots (small to challenge estimation)
true_angles = [-20 30]; % True DOAs (degrees)
d = 0.5;                % Element spacing
lambda = 1;
k = 2*pi/lambda;
SNR_dB_list = -20:5:20;
num_trials = 100;       % Monte Carlo trials

% Steering vector function
steering_vector = @(theta) exp(1j * k * d * (0:N-1).' * sind(theta));

% Preallocate error arrays
err_eig = zeros(size(SNR_dB_list));
err_svd = zeros(size(SNR_dB_list));
estimated_K_svd = zeros(size(SNR_dB_list));

for i = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(i);
    err_eig_trials = zeros(1, num_trials);
    err_svd_trials = zeros(1, num_trials);
    K_sum = 0;

    for trial = 1:num_trials
        %% === Generate highly correlated sources ===
        base_signal = randn(1, T) + 1j * randn(1, T);
        S = zeros(2, T);
        S(1, :) = base_signal;
        S(2, :) = base_signal + 0.00001 * (randn(1, T) + 1j * randn(1, T));  % Almost identical

        A = [steering_vector(true_angles(1)), steering_vector(true_angles(2))];
        noise = randn(N, T) + 1j * randn(N, T);
        signal_power = mean(abs(S(:)).^2);
        noise_power = signal_power / (10^(SNR_dB/10));
        X = A*S + sqrt(noise_power)*noise;

        %% === SVD-based ESPRIT with dynamic K ===
        [U, S_svd, ~] = svd(X, 'econ');
        sing_vals = diag(S_svd);
        threshold = 0.1 * max(sing_vals);
        K_svd = sum(sing_vals > threshold);
        K_svd = min(K_svd, N-1);
        K_sum = K_sum + K_svd;

        if K_svd < 1
            continue;
        end

        Us = U(:, 1:2);
        Phi_svd = (Us(1:end-1, :) \ Us(2:end, :));
        eigvals_svd = eig(Phi_svd);
        doas_svd = real(asind(angle(eigvals_svd) / (2*pi*d)));
        doas_svd = doas_svd(isfinite(doas_svd));
        if numel(doas_svd) >= 2
            doas_svd = sort(doas_svd(1:2));
            err_svd_trials(trial) = mean(abs(sort(true_angles(:)) - doas_svd(:)));
        end

        %% === EIG-based ESPRIT with fixed K = 2 ===
        Rxx = (X * X') / T;
        [V, D] = eig(Rxx);
        [~, idx] = sort(diag(D), 'descend');
        V = V(:, idx);
        Vs = V(:, 1:2);  % Fixed K = 2

        Phi_eig = (Vs(1:end-1, :) \ Vs(2:end, :));
        eigvals_eig = eig(Phi_eig);
        doas_eig = real(asind(angle(eigvals_eig) / (2*pi*d)));
        doas_eig = doas_eig(isfinite(doas_eig));
        if numel(doas_eig) >= 2
            doas_eig = sort(doas_eig(1:2));
            err_eig_trials(trial) = mean(abs(sort(true_angles(:)) - doas_eig(:)));
        end
    end

    % Average error over trials
    err_eig(i) = mean(err_eig_trials(err_eig_trials > 0));
    err_svd(i) = mean(err_svd_trials(err_svd_trials > 0));
    estimated_K_svd(i) = K_sum / num_trials;  % Mean K estimated via SVD
end

%% === Plot DOA estimation error ===
figure;
plot(SNR_dB_list, err_eig, 'r-o', 'LineWidth', 1.5); hold on;
plot(SNR_dB_list, err_svd, 'b--s', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('Mean DOA Estimation Error (degrees)');
title(['ESPRIT with Monte Carlo Averaging (' num2str(num_trials) ' trials)']);
legend('EIG-based (fixed K=2)', 'SVD-based (dynamic K)', 'Location', 'northeast');
grid on;

%% === Plot average estimated K from SVD ===
figure;
plot(SNR_dB_list, estimated_K_svd, 'k-x', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('Average Estimated K (via SVD)');
title('Rank Estimation from SVD');
grid on;
