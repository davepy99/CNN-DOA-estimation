clc; clear; close all;

%% Parameters
M = 8;
N = 200;
theta_true = 30;  % True DOA in degrees
lambda = 1;
d = lambda / 2;
k = 2*pi/lambda;
fft_len = 512;
trials = 100;
SNR_dB = [-10, 0, 10, 20, 30];
m = (0:M-1).';  % Antenna indices
theta_grid = linspace(-90, 90, fft_len);

%% Result storage
rmse_fft = zeros(size(SNR_dB));
rmse_music = zeros(size(SNR_dB));
rmse_esprit = zeros(size(SNR_dB));
rmse_mars = zeros(size(SNR_dB));

time_fft = zeros(size(SNR_dB));
time_music = zeros(size(SNR_dB));
time_esprit = zeros(size(SNR_dB));
time_mars = zeros(size(SNR_dB));

%% Loop over SNR
for snr_idx = 1:length(SNR_dB)
    snr = 10^(SNR_dB(snr_idx)/10);

    theta_fft_all = zeros(1, trials);
    theta_music_all = zeros(1, trials);
    theta_esprit_all = zeros(1, trials);
    theta_mars_all = zeros(1, trials);

    t_fft_all = zeros(1, trials);
    t_music_all = zeros(1, trials);
    t_esprit_all = zeros(1, trials);
    t_mars_all = zeros(1, trials);

    for t = 1:trials
        %% Signal simulation
        theta_rad = deg2rad(theta_true);
        a = exp(1j * k * d * m * sin(theta_rad));
        time_vec = 0:N-1;
        signal = exp(1j * 2 * pi * 0.05 * time_vec);
        X = a * signal;

        % Add noise
        noise = (randn(M, N) + 1j * randn(M, N)) / sqrt(2);
        X_noisy = X + noise * sqrt(norm(X,'fro')^2 / (norm(noise,'fro')^2 * snr));
        R = (1/N) * (X_noisy * X_noisy');

        %% FFT DOA (your version)
        tic;
        snapshot = X_noisy(:,1);
        fft_out = fftshift(fft(snapshot, fft_len));
        power_spectrum = abs(fft_out).^2;
        theta_axis = asind(linspace(-0.5, 0.5, fft_len) * lambda / d);
        [~, idx_fft] = max(power_spectrum);
        theta_fft_all(t) = theta_axis(idx_fft);
        t_fft_all(t) = toc;

        %% MUSIC
        tic;
        [V, D] = eig(R);
        [~, idx] = sort(diag(D), 'descend');
        En = V(:, idx(2:end));
        P_music = zeros(size(theta_grid));
        for i = 1:length(theta_grid)
            a_theta = exp(1j * k * d * m * sind(theta_grid(i)));
            P_music(i) = 1 / (a_theta' * (En * En') * a_theta);
        end
        [~, idx_music] = max(P_music);
        theta_music_all(t) = theta_grid(idx_music);
        t_music_all(t) = toc;

        %% ESPRIT
        tic;
        [Us, ~, ~] = svd(R);
        Us_sig = Us(:, 1);
        Us1 = Us_sig(1:end-1);
        Us2 = Us_sig(2:end);
        Phi = pinv(Us1) * Us2;
        eig_vals = eig(Phi);
        theta_esprit_all(t) = asind(angle(eig_vals(1)) / (k * d));
        t_esprit_all(t) = toc;

        %% MARS-like Estimation
        tic;
        cost_mars = zeros(size(theta_grid));
        for i = 1:length(theta_grid)
            a_theta = exp(1j * k * d * m * sind(theta_grid(i)));
            P = eye(M) - (a_theta * a_theta') / (a_theta' * a_theta);  % Orthogonal projection
            cost_mars(i) = real(trace(P * R));  % ML cost
        end
        [~, idx_mars] = min(cost_mars);  % ML: minimize cost
        theta_mars_all(t) = theta_grid(idx_mars);
        t_mars_all(t) = toc;
    end

    %% RMSE
    rmse_fft(snr_idx) = sqrt(mean((theta_fft_all - theta_true).^2));
    rmse_music(snr_idx) = sqrt(mean((theta_music_all - theta_true).^2));
    rmse_esprit(snr_idx) = sqrt(mean((theta_esprit_all - theta_true).^2));
    rmse_mars(snr_idx) = sqrt(mean((theta_mars_all - theta_true).^2));

    %% Time (ms)
    time_fft(snr_idx) = mean(t_fft_all) * 1000;
    time_music(snr_idx) = mean(t_music_all) * 1000;
    time_esprit(snr_idx) = mean(t_esprit_all) * 1000;
    time_mars(snr_idx) = mean(t_mars_all) * 1000;
end

%% Plot RMSE
figure;
plot(SNR_dB, rmse_fft, '-o', 'LineWidth', 2); hold on;
plot(SNR_dB, rmse_music, '-s', 'LineWidth', 2);
plot(SNR_dB, rmse_esprit, '-^', 'LineWidth', 2);
plot(SNR_dB, rmse_mars, '-x', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('RMSE (degrees)');
legend('FFT', 'MUSIC', 'ESPRIT', 'MARS-like');
title('DOA Estimation Accuracy (RMSE)');
grid on;



%% Plot RMSE w/o FFT
figure;
% plot(SNR_dB, rmse_fft, '-o', 'LineWidth', 2); 
plot(SNR_dB, rmse_music, '-s', 'LineWidth', 2);hold on;
plot(SNR_dB, rmse_esprit, '-^', 'LineWidth', 2);
plot(SNR_dB, rmse_mars, '-x', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('RMSE (degrees)');
legend('MUSIC', 'ESPRIT', 'MARS-like');
title('DOA Estimation Accuracy (RMSE)');
grid on;

%% Plot Computation Time
figure;
plot(SNR_dB, time_fft, '-o', 'LineWidth', 2); hold on;
plot(SNR_dB, time_music, '-s', 'LineWidth', 2);
plot(SNR_dB, time_esprit, '-^', 'LineWidth', 2);
plot(SNR_dB, time_mars, '-x', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Average Computation Time (ms)');
legend('FFT', 'MUSIC', 'ESPRIT', 'MARS-like');
title('Computation Time Comparison');
grid on;
