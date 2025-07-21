% Generate DOA test dataset for CNN evaluation

% Parameters (same structure as training)
fc = 1e9;
c = physconst('LightSpeed');
lambda = c / fc;
M = 8;
Nsamples = 1000;
Nexamples = 2000;                     % Smaller test set
SNR_dB = 10;                          % Match or vary from training
numSourcesRange = [1 3];
angleGrid = -90:2:90;                % Same resolution as training (2°)

% Array
array = phased.ULA('NumElements', M, 'ElementSpacing', lambda/2);
arrayPos = getElementPosition(array)/lambda;

% Data storage
covSet = zeros(M, M, 2*Nexamples);
labelSet = zeros(length(angleGrid), Nexamples);

for n = 1:Nexamples
    K = randi(numSourcesRange);
    doas = sort(rand(K,1)*180 - 90);  % Random DOAs

    A = steervec(arrayPos, doas.');
    s = randn(K, Nsamples);

    noise = (randn(M,Nsamples) + 1i*randn(M,Nsamples))/sqrt(2);
    SNR = 10^(SNR_dB/10);
    noisePower = norm(A*s,'fro')^2 / (SNR * norm(noise,'fro')^2);
    x = A*s + sqrt(noisePower)*noise;

    R = (x * x') / Nsamples;
    covSet(:,:,n) = real(R);
    covSet(:,:,n+Nexamples) = imag(R);

    for doa = doas'
        [~, idx] = min(abs(angleGrid - doa));
        labelSet(idx, n) = 1;
    end
end

% Save to new file
save('DOA_Test_Data_2deg.mat', 'covSet', 'labelSet', 'angleGrid', '-v7.3');
