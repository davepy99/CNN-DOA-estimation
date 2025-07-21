% Parameters
N = 8;             % Number of sensors
K = 3;             % Number of sources
snapshots = 200;   % Number of snapshots
SNR_dB = 10;       % Signal-to-noise ratio

% Array configuration
d = 0.5;           % Element spacing (in wavelength units)
angles = [-20, 0, 30];  % DOA of sources in degrees
lambda = 1;        % Wavelength
k = 2*pi/lambda;

% Steering matrix
A = zeros(N, K);
for i = 1:K
    A(:, i) = exp(1j * k * d * (0:N-1).' * sind(angles(i)));
end

% Signals and noise
S = randn(K, snapshots) + 1j*randn(K, snapshots);        % source signals
noise = (randn(N, snapshots) + 1j*randn(N, snapshots));  % noise
noise_power = norm(S, 'fro')^2 / (10^(SNR_dB/10)) / (N*snapshots);
X = A*S + sqrt(noise_power)*noise;                      % received signal

% Sample covariance matrix
Rxx = (X * X') / snapshots;

% Eigendecomposition
[Vecs, Vals] = eig(Rxx);
eigenvals = diag(Vals);
[eigenvals_sorted, idx] = sort(eigenvals, 'descend');
Vecs_sorted = Vecs(:, idx);

% Energy contribution visualization (eigenvalues)
figure;
stem(1:N, eigenvals_sorted, 'filled');
title('Eigenvalue Spectrum');
xlabel('Eigenvector Index');
ylabel('Eigenvalue (Energy)');
grid on;

% Optional: visualize cumulative energy
figure;
plot(cumsum(eigenvals_sorted)/sum(eigenvals_sorted), '-o');
title('Cumulative Energy Distribution');
xlabel('Number of Eigenvectors');
ylabel('Cumulative Energy');
grid on;
