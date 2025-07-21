clc; clear; close all;

%% Parameters
M = 8;                      % Number of array elements
d = 0.5;                    % Inter-element spacing (in wavelengths)
theta_true = 30;            % True DOA in degrees
lambda = 1;                 % Wavelength
k = 2*pi/lambda;            % Wave number
N_snapshots = 200;          % Number of snapshots
SNR_dB = [-10, 0, 10, 20, 30]; % SNR values to test
trials = 100;               % Monte Carlo trials per SNR

%% Result storage
rmse_fft = zeros(size(SNR_dB));
rmse_music = zeros(size(SNR_dB));
rmse_esprit = zeros(size(SNR_dB));

time_fft = zeros(size(SNR_dB));
time_music = zeros(size(SNR_dB));
time_esprit = zeros(size(SNR_dB));

%% DOA estimation over SNRs
for snr_idx = 1:length(SNR_dB)
    snr = 10^(SNR_dB(snr_idx)/10);
    theta_fft = zeros(1, trials);
    theta_music = zeros(1, trials);
    theta_esprit = zeros(1, trials);
    
    t_fft = zeros(1, trials);
    t_music = zeros(1, trials);
    t_esprit = zeros(1, trials);

    for t = 1:trials
        %% Signal Generation
        a = exp(1j*k*d*(0:M-1)'*sind(theta_true));
        signal = exp(1j*2*pi*rand(1, N_snapshots));
        x = a * signal;

        noise = (randn(M, N_snapshots) + 1j*randn(M, N_snapshots)) / sqrt(2);
        x_noisy = x + noise * sqrt(norm(x,'fro')^2 / (norm(noise,'fro')^2 * snr));
        R = (x_noisy * x_noisy') / N_snapshots;

        %% FFT DOA
        tic;
        fft_grid = -90:0.5:90;
        fft_resp = zeros(size(fft_grid));
        for i = 1:length(fft_grid)
            a_fft = exp(1j*k*d*(0:M-1)'*sind(fft_grid(i)));
            fft_resp(i) = abs(a_fft' * mean(x_noisy, 2));
        end
        [~, idx] = max(fft_resp);
        theta_fft(t) = fft_grid(idx);
        t_fft(t) = toc;

        %% MUSIC DOA
        tic;
        [V, D] = eig(R);
        [~, idx_sorted] = sort(diag(D), 'descend');
        En = V(:, idx_sorted(2:end));
        P_music = zeros(size(fft_grid));
        for i = 1:length(fft_grid)
            a_music = exp(1j*k*d*(0:M-1)'*sind(fft_grid(i)));
            P_music(i) = 1 / (a_music' * (En * En') * a_music);
        end
        [~, idx] = max(P_music);
        theta_music(t) = fft_grid(idx);
        t_music(t) = toc;

        %% ESPRIT DOA
        tic;
        [Us, ~, ~] = svd(R);
        Us_sig = Us(:, 1);
        Us1 = Us_sig(1:end-1);
        Us2 = Us_sig(2:end);
        Phi = pinv(Us1) * Us2;
        eig_vals = eig(Phi);
        theta_esprit(t) = asind(angle(eig_vals(1)) / (k * d));
        t_esprit(t) = toc;
    end

    %% RMSE Calculation
    rmse_fft(snr_idx) = sqrt(mean((theta_fft - theta_true).^2));
    rmse_music(snr_idx) = sqrt(mean((theta_music - theta_true).^2));
    rmse_esprit(snr_idx) = sqrt(mean((theta_esprit - theta_true).^2));

    %% Average Computation Time (per trial)
    time_fft(snr_idx) = mean(t_fft);
    time_music(snr_idx) = mean(t_music);
    time_esprit(snr_idx) = mean(t_esprit);
end

%% Plot RMSE
figure;
plot(SNR_dB, rmse_fft, '-o', 'LineWidth', 2); hold on;
plot(SNR_dB, rmse_music, '-s', 'LineWidth', 2);
plot(SNR_dB, rmse_esprit, '-^', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('DOA RMSE (degrees)');
legend('FFT', 'MUSIC', 'ESPRIT');
title('DOA Estimation Accuracy (RMSE)');

%% Plot Computation Time
figure;
plot(SNR_dB, time_fft * 1000, '-o', 'LineWidth', 2); hold on;
plot(SNR_dB, time_music * 1000, '-s', 'LineWidth', 2);
plot(SNR_dB, time_esprit * 1000, '-^', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Average Time per Trial (ms)');
legend('FFT', 'MUSIC', 'ESPRIT');
title('Computation Time Comparison');
