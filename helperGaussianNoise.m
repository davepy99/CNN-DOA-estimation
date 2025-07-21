function noise = helperGaussianNoise(M, SNR_dB, numSources)
    SNR_linear = 10^(SNR_dB/10);
    signal_power = 1;
    noise_power = signal_power / SNR_linear;
    noise = sqrt(noise_power/2) * (randn(M,1) + 1j*randn(M,1));
    
end
