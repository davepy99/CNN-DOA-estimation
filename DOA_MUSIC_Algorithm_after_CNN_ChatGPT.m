% --- Array Parameters ---
M = 8;                          % Number of antennas
d = 0.5;                        % Inter-element spacing in wavelengths
lambda = 1;                    % Wavelength
angle_scan = -90:1:90;         % Scan range in degrees
num_sources = numSources(1);   % Number of sources for snapshot 1

% --- Steering Vector Function ---
steerVec = @(theta) exp(1j * 2 * pi * d * (0:M-1)' * sind(theta));  % M×1 steering vector

% --- Eigen-decomposition of Covariance Matrix ---
R1 = R{1};                     % Covariance matrix for snapshot 1
[Evec, Eval] = eig(R1);
[evals_sorted, idx] = sort(diag(Eval), 'ascend');
En = Evec(:, idx(1:end-num_sources));  % Noise subspace

% --- MUSIC Spectrum Computation ---
music_spectrum = zeros(size(angle_scan));

for k = 1:length(angle_scan)
    a = steerVec(angle_scan(k));  % M×1
    music_spectrum(k) = 1 / (a' * (En * En') * a);
end

% --- Normalize and Plot Spectrum ---
music_spectrum_dB = 10*log10(abs(music_spectrum));
music_spectrum_dB = music_spectrum_dB - max(music_spectrum_dB);  % normalize

figure;
plot(angle_scan, music_spectrum_dB, 'LineWidth', 2);
xlabel('Angle (°)');
ylabel('Spatial Spectrum (dB)');
title('MUSIC Spectrum for R{1}');
grid on;

% --- Find Peaks (DOAs) ---
[~, peak_indices] = maxk(music_spectrum, num_sources);
estimated_DOAs = angle_scan(peak_indices);
disp('Estimated DOAs (degrees):');
disp(sort(estimated_DOAs));
