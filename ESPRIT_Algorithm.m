clc; clear; close all;

% === Parameters ===
M = 8;                      % Number of sensors
K = 5;                      % Number of sources
N = 1000;                   % Snapshots
fc = 1e9;                   % Frequency (Hz)
c = 3e8;                    % Speed of light
lambda = c/fc;              % Wavelength
d = lambda/2;               % Element spacing
                % SNR in dB
trueAngles = [10, 20, 30, 40, 50];      % Degrees

figure;
hold on 
SNR_dB1 = [-10 0 10 20];
colors = lines(numel(SNR_dB1));

f= 0;
for SNR_dB = SNR_dB1
f = f+1;
% === Steering Matrix ===
angles_rad = deg2rad(trueAngles);
a = @(theta) exp(1j*2*pi*d*sin(theta)/lambda * (0:M-1).'); % steering vector
A = [a(angles_rad(1)), a(angles_rad(2)), a(angles_rad(3)), a(angles_rad(4)), a(angles_rad(5))];                  % M x K

% === Signal + Noise ===
S = (randn(K,N) + 1i*randn(K,N))/sqrt(2);                  % Random source signals
X = A * S;                                                 % Ideal array signal
sigma2 = 10^(-SNR_dB/10);
X = X + sqrt(sigma2/2)*(randn(M,N) + 1i*randn(M,N));       % Add complex Gaussian noise

% === Covariance and Subspace ===
R = (X*X')/N;
[U,~,~] = svd(X);                     % Singular Value Decomposition (SVD)
Us = U(:,1:K);                        % Signal subspace (M x K)

[R_V, R_D] = eig(R);
[R_eigenvals, idx] = sort(diag(R_D), 'descend');
R_V = R_V(:, idx);                    % Sort eigenvectors by descending eigenvalues
Vs = R_V(:, 1:K);                     % Signal subspace

% === Subarray Extraction ===
% You can use Vs or Us, one is using eig() the other is unsing svd()
Us1 = Us(1:end-1, :);                % First subarray
Us2 = Us(2:end, :);                  % Second subarray

% === Solve Rotation Matrix ===
Phi = pinv(Us1) * Us2;
eigvals = eig(Phi);
phi_hat = angle(eigvals);           % Phase shifts

% === Estimate DOAs ===
theta_hat = asin(phi_hat * lambda / (2*pi*d)) * 180/pi;

% === Display Results ===
disp('True DOAs (degrees):');
disp(sort(trueAngles));

disp('Estimated DOAs (degrees):');
disp(sort(real(theta_hat)));

% === Plot array response (optional) ===
scanAngles = -90:0.5:90;
a_scan = @(theta) exp(1j*2*pi*d*sind(theta)/lambda * (0:M-1).');
P = zeros(size(scanAngles));
for k = 1:length(scanAngles)
    av = a_scan(scanAngles(k));
    P(k) = 1 / (av' * (eye(M) - Us*Us') * av);
end
P_dB = 10*log10(abs(P)/max(abs(P)));





% === Subarray Extraction2 ===
% You can use Vs or Us, one is using eig() the other is unsing svd()
Us1_2 = Vs(1:end-1, :);                % First subarray
Us2_2 = Vs(2:end, :);                  % Second subarray

% === Solve Rotation Matrix ===
Phi_2 = pinv(Us1_2) * Us2_2;
eigvals_2 = eig(Phi_2);
phi_hat_2 = angle(eigvals_2);           % Phase shifts

% === Estimate DOAs ===
theta_hat_2 = asin(phi_hat_2 * lambda / (2*pi*d)) * 180/pi;

% === Display Results ===
disp('True DOAs Number 2 (degrees):');
disp(sort(trueAngles));

disp('Estimated DOAs Number 2(degrees):');
disp(sort(real(theta_hat_2)));

% === Plot array response (optional) ===
scanAngles = -90:0.5:90;
a_scan = @(theta) exp(1j*2*pi*d*sind(theta)/lambda * (0:M-1).');
P_2 = zeros(size(scanAngles));
for k = 1:length(scanAngles)
    av = a_scan(scanAngles(k));
    P_2(k) = 1 / (av' * (eye(M) - Vs*Vs') * av);
end
P_dB_2 = 10*log10(abs(P_2)/max(abs(P_2)));






plot(scanAngles, P_dB,  '--', 'Color', colors(f, :), 'LineWidth', 2, 'HandleVisibility', 'off'); 
plot(scanAngles, P_dB_2,  ':', 'Color', colors(f, :), 'LineWidth', 2, 'HandleVisibility', 'off'); 

xline(real(theta_hat), '--', 'Color', colors(f, :),  'LineWidth', 1.5, 'HandleVisibility', 'off');
h= xline(real(theta_hat_2), '--', 'Color', colors(f, :), 'LineWidth', 1.5);
set(h, 'DisplayName', ['Cycle ' num2str(f)]);
end


xline(trueAngles, 'black', 'LineWidth', 3, 'HandleVisibility', 'off');
xlabel('Angle (degrees)');
ylabel('Pseudo-Spectrum (dB)');
legend show;
title('ESPRIT-Based DOA Estimation');
grid on;
hold off;
