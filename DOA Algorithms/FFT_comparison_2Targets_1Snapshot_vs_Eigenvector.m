clc; clear; close all;

% Parameters
M = 8;
N = 200;
lambda = 1;
d = lambda / 2;
k = 2 * pi / lambda;
t = 0:N-1;
fft_len = 512;
fs = linspace(-0.5, 0.5, fft_len);
fs_clipped = min(max(fs * lambda / d, -1), 1);
theta_axis = asind(fs_clipped);

% Two source DOAs
theta1 = 20;
theta2 = -30;
theta1_rad = deg2rad(theta1);
theta2_rad = deg2rad(theta2);

% Steering vectors
m = (0:M-1).';
a1 = exp(1j * k * d * m * sin(theta1_rad));
a2 = exp(1j * k * d * m * sin(theta2_rad));

% Signals
s1 = exp(1j * 2 * pi * 0.05 * t);
s2 = exp(1j * 2 * pi * 0.12 * t);

% Combine signals
X_clean = a1 * s1 + a2 * s2;

% SNR sweep
snr_dBs = [0 10 20];
colors = lines(length(snr_dBs));
figure;

for i = 1:length(snr_dBs)
    snr_dB = snr_dBs(i);
    snr_linear = 10^(snr_dB / 10);
    noise_power = 1 / snr_linear;

    % Add noise
    noise = sqrt(noise_power/2) * (randn(M, N) + 1j * randn(M, N));
    X_noisy = X_clean + noise;

    % --- Algorithm 1: Dominant eigenvector FFT
    R = (1/N) * (X_noisy * X_noisy');
    [V, D] = eig(R);
    [~, idxs] = sort(real(diag(D)), 'descend');
    v1 = V(:, idxs(1));
    v2 = V(:, idxs(2));
    fft1 = fftshift(fft(v1, fft_len));
    fft2 = fftshift(fft(v2, fft_len));
    power1 = abs(fft1).^2 + abs(fft2).^2;  % combine both modes
    power1_dB = 10*log10(power1 + 1e-12);
    power1_dB = power1_dB - max(power1_dB);

    % --- Algorithm 2: Average snapshot
    snapshot_avg = mean(X_noisy, 2);
    fft_snapshot = fftshift(fft(snapshot_avg, fft_len));
    power2 = abs(fft_snapshot).^2;
    power2_dB = 10*log10(power2 + 1e-12);
    power2_dB = power2_dB - max(power2_dB);

    % Plot both
    plot(theta_axis, power1_dB, '-', 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', ['Eigenvec SNR ' num2str(snr_dB) ' dB']); hold on;
    % plot(theta_axis, power2_dB, '--', 'Color', colors(i,:), 'LineWidth', 1.5, ...
    %     'DisplayName', ['Snapshot SNR ' num2str(snr_dB) ' dB']);
end

% True DOAs
xline(theta1, ':r', 'LineWidth', 1.5, 'DisplayName', 'True DOA 1');
xline(theta2, ':b', 'LineWidth', 1.5, 'DisplayName', 'True DOA 2');

xlabel('Angle (degrees)');
ylabel('Normalized Power (dB)');
title('FFT-Based DOA Estimation for Two Sources');
legend('Location', 'best');
grid on;
xlim([-90 90]);
