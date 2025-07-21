clc; clear; close all;

%% Parameters
M = 16;                  % Number of antenna elements
d = 0.5;                 % Element spacing (in wavelengths)
lambda = 1;             % Wavelength
theta_true = 20;        % True DOA in degrees
k = 2 * pi / lambda;    % Wave number
N_snapshots = 200;      % Number of time samples

%% Generate incoming signal
m = (0:M-1).';                          % Sensor indices
theta_rad = deg2rad(theta_true);
a_theta = exp(1j * k * d * m * sin(theta_rad));  % Steering vector

signal = exp(1j * 2 * pi * 0.05 * (0:N_snapshots-1));  % Time-domain signal
X = a_theta * signal;                    % M x N received signal (ideal)

%% Add noise
SNR_dB = 10;
SNR = 10^(SNR_dB/10);
noise = (randn(M, N_snapshots) + 1j * randn(M, N_snapshots)) / sqrt(2);
X_noisy = X + noise * sqrt(norm(X,'fro')^2 / (norm(noise,'fro')^2 * SNR));

%% Snapshot average
snapshot = mean(X_noisy, 2);   % M x 1

%% FFT across array elements
fft_len = 512;
spatial_fft = fftshift(fft(snapshot, fft_len));
power_spectrum = abs(spatial_fft).^2;

%% Angle mapping (from spatial frequency)
fs = linspace(-0.5, 0.5, fft_len);
theta_axis = asind(fs * lambda / d);

%% Plot FFT Spatial Spectrum
figure;
plot(theta_axis, 10*log10(power_spectrum / max(power_spectrum)), 'LineWidth', 2);
xlabel('Angle (degrees)');
ylabel('Normalized Power (dB)');
title('FFT-Based DOA Estimation');
grid on;
xlim([-90, 90]);
hold on;
xline(theta_true, '--r', 'True DOA');
legend('FFT Spectrum', 'True DOA');
