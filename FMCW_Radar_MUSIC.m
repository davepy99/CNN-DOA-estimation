clc; clear; close all;

%% === RADAR & TARGET PARAMETERS ===
c = 3e8;             % Speed of light (m/s)
fc = 77e9;           % Carrier frequency (77 GHz)
lambda = c / fc;

% FMCW chirp parameters
B = 50e6;            % Bandwidth (50 MHz)
T_chirp = 20e-6;     % Chirp duration (20 us)
slope = B / T_chirp; % Chirp slope (Hz/s)

% Sampling
fs = 5e6;                         % Lowered sampling rate for memory safety
N_fast = round(T_chirp * fs);    % Fast time samples per chirp
N_slow = 64;                     % Number of chirps (slow time)

% Array configuration
M = 8;                           % Number of antenna elements
d = lambda / 2;

% Targets: [range (m), velocity (m/s), angle (deg)]
targets = [
    30, 20, 20;
    40, -10, 50
];
K = size(targets,1);

%% === TIME AXES ===
t_fast = (0:N_fast-1) / fs;       % fast time within chirp
t_slow = (0:N_slow-1) * T_chirp;  % slow time (chirp index)

%% === GENERATE SIGNALS AT EACH ANTENNA ===
X = zeros(N_fast, N_slow, M);     % radar cube: fast × slow × antenna

for m = 1:M
    for k = 1:K
        R = targets(k,1); V = targets(k,2); theta = targets(k,3);
        phi_m = 2 * pi * d * (m-1) * sind(theta) / lambda;

        for n = 1:N_slow
            Rn = R + V * t_slow(n);
            tau_n = 2 * Rn / c;
            fb = slope * tau_n;
            fd = 2 * V / lambda;

            beat = exp(1j*2*pi*(fb * t_fast + fd * t_slow(n)) + 1j * phi_m);
            X(:,n,m) = X(:,n,m) + beat(:);  % Ensure beat is column
        end
    end
end

% Add noise
SNR_dB = 20;
noise_power = 10^(-SNR_dB/10);
X = X + sqrt(noise_power/2)*(randn(size(X)) + 1i*randn(size(X)));

%% === RANGE & DOPPLER FFT ===
X_r = fft(X, [], 1);                  % Range FFT
X_r = X_r(1:N_fast/2,:,:);            % Keep positive frequencies
N_range = size(X_r,1);

X_rd = fftshift(fft(X_r, [], 2),2);   % Doppler FFT
N_doppler = size(X_rd,2);

%% === AUTO BIN SELECTION BASED ON ENERGY ===
energyMap = squeeze(sum(abs(X_rd).^2, 3)); % Sum over antennas
[maxVal, idx] = max(energyMap(:));
[range_bin, doppler_bin] = ind2sub(size(energyMap), idx);
fprintf('Selected bin: Range=%d, Doppler=%d, Energy=%.2f\n', range_bin, doppler_bin, maxVal);

% Plot range-Doppler map
figure;
imagesc(-N_doppler/2:N_doppler/2-1, 1:N_range, 10*log10(energyMap));
xlabel('Doppler bin'); ylabel('Range bin'); title('Range-Doppler Map (dB)');
colorbar;

%% === NEW MUSIC SECTION: MULTIPLE SNAPSHOTS FROM CHIRPS (SLOW-TIME) ===
% Extract range bin across all slow-time samples (chirps)
X_snapshot = squeeze(X_r(range_bin, :, :));  % [slow-time x antennas] = [N_slow x M]
X_snapshot = X_snapshot.';                  % Transpose → [M x N_slow]

% Covariance matrix from slow-time samples
R = (X_snapshot * X_snapshot') / N_slow;

% MUSIC subspace
[U,~,~] = svd(R);
K_est = 2; % Known number of sources
Un = U(:, K_est+1:end);

% MUSIC spatial spectrum
angles = -90:0.5:90;
Pmusic = zeros(size(angles));
for i = 1:length(angles)
    a = exp(1j * 2 * pi * d * sind(angles(i)) / lambda * (0:M-1)).';
    Pmusic(i) = 1 / (a' * (Un * Un') * a);
end
Pmusic_dB = 10 * log10(abs(Pmusic) / max(abs(Pmusic)));

% Plot MUSIC spectrum
figure;
plot(angles, Pmusic_dB, 'b', 'LineWidth', 2); grid on;
xlabel('Angle (degrees)');
ylabel('Pseudo-spectrum (dB)');
title('MUSIC DOA Estimation Using Slow-Time Snapshots');
xlim([-90 90]);

drawnow;                % Force graphics update
pause(0.1);            % Let GUI thread catch up
