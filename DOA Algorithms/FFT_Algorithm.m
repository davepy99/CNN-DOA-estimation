clc; clear; close all;

%% Parameters
M = 8;                      % Number of antenna elements
N = 200;                    % Number of time snapshots
theta_true = 30;            % True DOA in degrees
lambda = 1;                 % Wavelength
d = lambda / 2;             % Inter-element spacing
k = 2 * pi / lambda;        % Wave number
fft_len = 512;              % FFT length
trials = 100;               % Monte Carlo trials per SNR
SNR_dB = [-10, 0, 10, 20, 30];

%% Frequency to angle axis for FFT
fs = linspace(-0.5, 0.5, fft_len);
theta_axis = asind(fs * lambda / d);

%% Result storage
rmse_fft = zeros(size(SNR_dB));
rmse_music = zeros(size(SNR_dB));
rmse_esprit = zeros(size(SNR_dB));

time_fft = zeros(size(SNR_dB));
time_music = zeros(size(SNR_dB));
time_esprit = zeros(size(SNR_dB));

%% Loop over SNRs
for snr_idx = 1:length(SNR_dB)
    snr = 10^(SNR_dB(snr_idx)/10);
    theta_fft_all = zeros(1, trials);
    theta_music_all = zeros(1, trials);
    theta_esprit_all = zeros(1, trials);

    t_fft_all = zeros(1, trials);
    t_music_all = zeros(1, trials);
    t_esprit_all = zeros(1, trials);

    for t = 1:trials
        %% Signal generation
        m = (0:M-1).';  % Antenna indices
        theta_rad = deg2rad(theta_true);
        a = exp(1j * k * d * m * sin(theta_rad));  % Steering vector
        time_vec = 0:N-1;
        signal = exp(1j * 2 * pi * 0.05 * time_vec);  % Narrowband signal
        X = a * signal;  % M x N ideal signal

        % Add noise
        noise_power = 1 / snr;
        noise = sqrt(noise_power/2) * (randn(M, N) + 1j * randn(M, N));
        X_noisy = X + noise;

        %% FFT-Based DOA Estimation (your method)
        tic;
        snapshot = X_noisy(:,1);  % One time snapshot
        fft_out = fftshift(fft(snapshot, fft_len));
        power_spectrum = abs(fft_out).^2;
        [~, idx_fft] = max(power_spectrum);
        theta_fft_all(t) = theta_axis(idx_fft);
        t_fft_all(t) = toc;

        %% MUSIC-Based DOA Estimation
        tic;
        R = (1/N) * (X_noisy * X_noisy');
        [V, D] = eig(R);
        [~, idx] = sort(diag(D), 'descend');
        En = V(:, idx(2:end));  % Noise subspace
        P_music = zeros(size(theta_axis));
        for i = 1:length(theta_axis)
            a_theta = exp(1j * k * d * m * sind(theta_axis(i)));
            P_music(i) = 1 / (a_theta' * (En * En') * a_theta);
        end
        [~, idx_music] = max(P_music);
        theta_music_all(t) = theta_axis(idx_music);
        t_music_all(t) = toc;

        %% ESPRIT-Based DOA Estimation
        tic;
        [Us, ~, ~] = svd(R);
        Us_sig = Us(:, 1);
        Us1 = Us_sig(1:end-1);
        Us2 = Us_sig(2:end);
        Phi = pinv(Us1) * Us2;
        eig_vals = eig(Phi);
        theta_esprit_all(t) = asind(angle(eig_vals(1)) / (k * d));
        t_esprit_all(t) = toc;
    end

    %% RMSE
    rmse_fft(snr_idx) = sqrt(mean((theta_fft_all - theta_true).^2));
    rmse_music(snr_idx) = sqrt(mean((theta_music_all - theta_true).^2));
    rmse_esprit(snr_idx) = sqrt(mean((theta_esprit_all - theta_true).^2));

    %% Timing
    time_fft(snr_idx) = mean(t_fft_all) * 1000;     % ms
    time_music(snr_idx) = mean(t_music_all) * 1000;
    time_esprit(snr_idx) = mean(t_esprit_all) * 1000;
end

%% Plot RMSE
figure;
plot(SNR_dB, rmse_fft, '-o', 'LineWidth', 2); hold on;
plot(SNR_dB, rmse_music, '-s', 'LineWidth', 2);
plot(SNR_dB, rmse_esprit, '-^', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('RMSE (degrees)');
legend('FFT', 'MUSIC', 'ESPRIT');
title('DOA Estimation Accuracy');
grid on;

%% Plot Timing
figure;
plot(SNR_dB, time_fft, '-o', 'LineWidth', 2); hold on;
plot(SNR_dB, time_music, '-s', 'LineWidth', 2);
plot(SNR_dB, time_esprit, '-^', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Average Computation Time per Trial (ms)');
legend('FFT', 'MUSIC', 'ESPRIT');
title('Computation Time Comparison');
grid on;
