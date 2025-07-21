clc; clear; close all;

%% === PARAMETERS ===
M = 8;                        % Number of antenna elements
K = 2;                        % Number of sources
N = 50;                    % Number of snapshots
fc = 1e9;                    % Carrier frequency (Hz)
c = 3e8;                     % Speed of light (m/s)
lambda = c / fc;             % Wavelength
d = lambda / 2;              % Element spacing (half-wavelength)
SNR_dB = 20;                 % SNR in dB
trueAngles = [0 30];       % True DOAs in degrees
maxNumCompThreads(1);
%% === SIGNAL MODEL ===
% Generate steering matrix A
a = @(theta) exp(1j * 2 * pi * d * sind(theta) / lambda * (0:M-1).'); % Mx1 steering vector
A = [a(trueAngles(1)), a(trueAngles(2))];  % MxK

% Generate random source signals (complex baseband)
S = (randn(K,N) + 1i * randn(K,N)) / sqrt(2);  % KxN

% Generate received signal X = A*S + noise
X = A * S;
sigma2 = 10^(-SNR_dB / 10);
X = X + sqrt(sigma2/2) * (randn(M,N) + 1i*randn(M,N));  % AWGN

%% === STEP 1: Compute Covariance Matrix ===
tic
R = (X * X') / N;

%% === STEP 2: Eigen-Decomposition ===
[U, D] = eig(R);
[eigvals, idx] = sort(diag(D), 'descend');
U = U(:, idx);  % Sort eigenvectors accordingly

% Separate signal and noise subspaces
Us = U(:, 1:K);          % Signal subspace
Un = U(:, K+1:end);      % Noise subspace
toc
%% === STEP 3: MUSIC Spatial Spectrum ===
tic
scanAngles = -5:0.5:5;
% scanAngles = -90:0.5:90;
Pmusic = zeros(size(scanAngles));

for i = 1:length(scanAngles)
    steeringVec = a(scanAngles(i));
    Pmusic(i) = 1 / (steeringVec' * (Un * Un') * steeringVec);
end
toc
% Convert to dB scale
Pmusic_dB = 10 * log10(abs(Pmusic) / max(abs(Pmusic)));

%% === STEP 4: Plot MUSIC Spectrum ===
figure;
plot(scanAngles, Pmusic_dB, 'b', 'LineWidth', 2); hold on;
xline(trueAngles(1), 'k--', 'LineWidth', 1.5);
xline(trueAngles(2), 'k--', 'LineWidth', 1.5);
xlabel('Angle (degrees)');
ylabel('Pseudo-Spectrum (dB)');
title('MUSIC DOA Estimation');
legend('MUSIC Spectrum', 'True DOAs');
grid on;
xlim([-90, 90]);
