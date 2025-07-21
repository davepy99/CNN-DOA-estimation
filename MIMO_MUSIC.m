% === Configuration ===
fc = 1e9;                                % 1 GHz carrier frequency
c = physconst('LightSpeed');
lambda = c / fc;

% === Tx and Rx Arrays ===
txArray = phased.ULA('NumElements', 4, 'ElementSpacing', lambda/2);
rxArray = phased.ULA('NumElements', 4, 'ElementSpacing', lambda/2);

% === Virtual Array Formation ===
txPos = txArray.getElementPosition();
rxPos = rxArray.getElementPosition();

virtualPos = [];
for m = 1:size(txPos, 2)
    for n = 1:size(rxPos, 2)
        virtualPos(:, end+1) = txPos(:, m) + rxPos(:, n);  % virtual element = tx + rx
    end
end

% === Target DOAs ===
doas = [30, -20];              % Two incoming targets
numTargets = numel(doas);
Nsamples = 1000;               % Time samples

% === Orthogonal Waveforms from Tx ===
waveforms = eye(4);            % Orthogonal waveforms (identity matrix for simplicity)

% === Signal Generation ===
sig = randn(Nsamples, numTargets);   % Baseband signals per target

% === Create Received Signals at Virtual Elements ===
rxCombined = zeros(size(virtualPos, 2), Nsamples);  % virtual array data

% For each target
for k = 1:numTargets
    angle = doas(k);
    
    % Steering vectors for Tx and Rx arrays
    atx = steervec(txPos/lambda, angle);  % [M x 1]
    arx = steervec(rxPos/lambda, angle);  % [N x 1]
    
    % Outer product: virtual array response [M*N x 1]
    avirtual = kron(arx, atx);  % Vectorized response of virtual array
    
    % Received signal from this target
    targetSignal = sig(:, k).';  % Row vector [1 x Nsamples]
    
    % Add target contribution to total signal
    rxCombined = rxCombined + avirtual * targetSignal;
end

% === Add Noise (optional) ===
rxCombined = rxCombined + 0.01*randn(size(rxCombined));

% === Estimate DOAs with MUSIC ===
% Compute covariance matrix
R = (rxCombined * rxCombined') / Nsamples;

% Eigendecomposition
[Evec, Eval] = eig(R);
[~, idx] = sort(diag(Eval), 'descend');
En = Evec(:, numTargets+1:end);  % Noise subspace

% MUSIC spectrum
scanAngles = -90:0.5:90;
spectrum = zeros(size(scanAngles));

for k = 1:length(scanAngles)
    a = steervec(virtualPos/lambda, scanAngles(k));
    spectrum(k) = 1 / abs(a' * En * En' * a);
end

% === Normalize and Plot ===
spectrum = real(spectrum);
spectrum_dB = 10 * log10(spectrum / max(spectrum));  % normalize to 0 dB peak

figure;
plot(scanAngles, spectrum_dB, 'LineWidth', 1.5);
xlabel('Angle (degrees)');
ylabel('Spatial Spectrum (dB)');
title('MIMO MUSIC Spectrum using Virtual Array');
grid on;
