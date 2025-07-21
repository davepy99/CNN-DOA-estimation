% Parameters
fc = 1e9;                         % Carrier frequency (Hz)
c = physconst('LightSpeed');     % Speed of light
lambda = c/fc;
M = 8;                            % Number of array elements
Nsamples = 1000;                  % Snapshots per sample
Nexamples = 20000;                % Increased dataset size
SNR_dB = 10;                      % SNR in dB
numSourcesRange = [1 3];          % Number of sources per sample
angleGrid = -90:2:90;             % 2-degree resolution (91 classes)

% Array
array = phased.ULA('NumElements', M, 'ElementSpacing', lambda/2);
arrayPos = getElementPosition(array)/lambda;

% Dataset containers
covSet = zeros(M, M, 2*Nexamples);           % Real & imag split
labelSet = zeros(length(angleGrid), Nexamples);  % Binary DOA label vector

for n = 1:Nexamples
    K = randi(numSourcesRange);
    doas = sort(rand(K,1)*180 - 90);  % Random DOAs between -90° and 90°
    
    % Steering matrix & signal
    A = steervec(arrayPos, doas.');
    s = randn(K, Nsamples);
    
    noise = (randn(M,Nsamples) + 1i*randn(M,Nsamples))/sqrt(2);
    SNR = 10^(SNR_dB/10);
    noisePower = norm(A*s,'fro')^2 / (SNR * norm(noise,'fro')^2);
    x = A*s + sqrt(noisePower)*noise;
    
    % Covariance
    R = (x * x') / Nsamples;
    covSet(:,:,n)       = real(R);
    covSet(:,:,n+Nexamples) = imag(R);
    
    % Label (binary vector over 2° bins)
    for doa = doas'
        [~, idx] = min(abs(angleGrid - doa));
        labelSet(idx, n) = 1;
    end
end

% Save data
save('DOA_Training_Data_2deg.mat', 'covSet', 'labelSet', 'angleGrid', '-v7.3');
