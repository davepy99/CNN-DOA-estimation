clc; clear; close all;

% Parameters
M = 8;
N = 200;
theta_deg = 25;
lambda = 1;
d = lambda / 2;
k = 2 * pi / lambda;
t = 0:N-1;

% Signal and steering
signal = exp(1j * 2 * pi * 0.05 * t);
m = (0:M-1).';
theta_rad = deg2rad(theta_deg);
a_theta = exp(1j * k * d * m * sin(theta_rad));

% FFT setup
fft_len = 512;
fs = linspace(-0.5, 0.5, fft_len);
fs_clipped = min(max(fs * lambda / d, -1), 1);  % avoid NaNs in asind
theta_axis = asind(fs_clipped);

% SNRs
snr_dBs = [-20 -10 0 10 20 30];
colors = lines(length(snr_dBs));
figure;

for i = 1:length(snr_dBs)
    snr_dB = snr_dBs(i);
    snr_linear = 10^(snr_dB / 10);
    noise_power = 1 / snr_linear;

    % Simulate signal with noise
    X_clean = a_theta * signal;
    noise = sqrt(noise_power / 2) * (randn(M, N) + 1j * randn(M, N));
    X_noisy = X_clean + noise;

    % --- Algorithm 1: Eigenvector FFT
    R = (1/N) * (X_noisy * X_noisy');
    [V, D] = eig(R);
    [~, idx] = max(real(diag(D)));
    dominant = V(:, idx);
    fft1 = fftshift(fft(dominant, fft_len));
    power1 = 10 * log10(abs(fft1).^2 + 1e-12);  % add epsilon to avoid -Inf
    power1 = power1 - max(power1);

    % --- Algorithm 2: Snapshot average FFT
    avg_snapshot = X_noisy(:,1);
    fft2 = fftshift(fft(avg_snapshot, fft_len));
    power2 = 10 * log10(abs(fft2).^2 + 1e-12);
    power2 = power2 - max(power2);

    % Plot both
    plot(theta_axis, power1, '-', 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', ['Alg1 SNR ' num2str(snr_dB) ' dB']); hold on;
    plot(theta_axis, power2, '--', 'LineWidth', 1.5, 'Color', colors(i,:), ...
        'DisplayName', ['Alg2 SNR ' num2str(snr_dB) ' dB']);
end

xlabel('Angle (degrees)');
ylabel('Normalized Power (dB)');
title('FFT-Based DOA vs SNR (Two Algorithms)');
legend('Location', 'best');
grid on;
xlim([-90 90]);
