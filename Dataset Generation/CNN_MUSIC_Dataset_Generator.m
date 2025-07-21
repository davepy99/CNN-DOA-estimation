% Parameters
fc = 1e9;                      % Carrier frequency (Hz)
c = physconst('LightSpeed');  % Speed of light
lambda = c/fc;
M = 8;                        % Number of array elements
Nsamples = 1000;              % Snapshot length
Nexamples = 5000;             % Total training examples
SNR_dB = 10;                  % SNR in dB
numSourcesRange = [1 3];      % Random number of sources per sample
angleGrid = -90:1:90;         % Label binning resolution (for classification)

% Array setup
array = phased.ULA('NumElements', M, 'ElementSpacing', lambda/2);
arrayPos = getElementPosition(array)/lambda;

% Initialize dataset
covSet = zeros(M, M, Nexamples);  % complex covariance matrices
labelSet = zeros(length(angleGrid), Nexamples);  % multi-label output (1s where sources exist)

% Main loop
for n = 1:Nexamples
    % Number of sources
    K = randi(numSourcesRange);
    
    % True DOAs
    doas = sort(rand(K,1)*180 - 90);  % Uniformly in [-90, 90]
    
    % Signal generation
    A = steervec(arrayPos, doas.');
    s = randn(K, Nsamples);  % baseband signals (can be improved)
    
    % Noise
    noise = (randn(M,Nsamples) + 1i*randn(M,Nsamples))/sqrt(2);
    SNR = 10^(SNR_dB/10);
    noisePower = norm(A*s,'fro')^2 / (SNR * norm(noise,'fro')^2);
    x = A*s + sqrt(noisePower)*noise;
    
    % Sample covariance matrix
    R = (x * x') / Nsamples;
    covSet(:,:,n) = R;

    % Label: multi-label angle classification (binary vector over angle grid)
    for doa = doas'
        [~, idx] = min(abs(angleGrid - doa));
        labelSet(idx, n) = 1;
    end
end

% Optionally separate into real and imag parts for CNN input
covSet_real = real(covSet);
covSet_imag = imag(covSet);
X = cat(3, covSet_real, covSet_imag);  % [M x M x 2 x Nexamples] if needed

% Save dataset
save('DOA_Training_Data.mat', 'X', 'labelSet', 'angleGrid', '-v7.3');
