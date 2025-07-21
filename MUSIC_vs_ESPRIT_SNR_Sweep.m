clc; clear; close all;

% === Parameters ===
fc = 1e9;                            % Carrier frequency (Hz)
c = 3e8;                             % Speed of light (m/s)
lambda = c/fc;                       % Wavelength
d = lambda/2;                        % Element spacing
M = 8;                               % Number of antenna elements
N = 1000;                            % Number of snapshots
trueAngles = [20 50];               % True DOAs in degrees
K = length(trueAngles);             % Number of sources
numTrials = 100;                    % Monte Carlo runs per SNR
snrVals = -10:2:20;                 % SNR values to sweep (dB)
RMSE_MUSIC = zeros(size(snrVals));
RMSE_ESPRIT = zeros(size(snrVals));

% === ULA and steering vector config ===
ULA = phased.ULA('NumElements', M, 'ElementSpacing', d);
steeringVec = phased.SteeringVector('SensorArray', ULA, ...
    'IncludeElementResponse', false, 'PropagationSpeed', c);

% === SNR Sweep ===
for s = 1:length(snrVals)
    SNR_dB = snrVals(s);
    noisePower = 10^(-SNR_dB/10);
    err_music = zeros(numTrials,1);
    err_esprit = zeros(numTrials,1);

    for t = 1:numTrials
        % Simulate signal
        A = steeringVec(fc, trueAngles);
        S = (randn(K,N) + 1i*randn(K,N))/sqrt(2);    % Random complex signals
        X = A*S;

        % Add noise
        noise = sqrt(noisePower/2)*(randn(M,N) + 1i*randn(M,N));
        Xn = X + noise;

        % MUSIC Estimation
        R = Xn*Xn'/N;
        [U,~,~] = svd(R);
        Un = U(:,K+1:end);  % Noise subspace
        scanAngles = -90:0.5:90;
        Pmusic = zeros(size(scanAngles));
        for ii = 1:length(scanAngles)
            a = steeringVec(fc, scanAngles(ii));
            Pmusic(ii) = 1 / (a' * (Un * Un') * a);
        end
        [~, locs] = findpeaks(real(Pmusic), 'SortStr', 'descend', 'NPeaks', K);
        estAngles_music = sort(scanAngles(locs(1:K)));
        err_music(t) = sqrt(mean((sort(trueAngles) - estAngles_music).^2));

        % ESPRIT Estimation (robust)
        try
            [Ue,~,~] = svd(R);
            Us = Ue(:,1:K);
            Us1 = Us(1:end-1,:);
            Us2 = Us(2:end,:);
            Phi = pinv(Us1)*Us2;
            eigvals = eig(Phi);
            doa_esprit_all = asin(angle(eigvals)*lambda/(2*pi*d)) * 180/pi;
            doa_esprit_all = real(doa_esprit_all(~isnan(doa_esprit_all) & isreal(doa_esprit_all)));

            % Choose closest K DOAs to the true ones
            if length(doa_esprit_all) >= K
                matched = zeros(1,K);
                for k = 1:K
                    [~, idx] = min(abs(doa_esprit_all - trueAngles(k)));
                    matched(k) = doa_esprit_all(idx);
                    doa_esprit_all(idx) = Inf;  % prevent reuse
                end
                doa_esprit = sort(matched);
                err_esprit(t) = sqrt(mean((sort(trueAngles) - doa_esprit).^2));
            else
                err_esprit(t) = NaN;
            end
        catch
            err_esprit(t) = NaN;
        end
    end

    % Store average RMSE (excluding failed ESPRIT trials)
    RMSE_MUSIC(s) = mean(err_music);
    RMSE_ESPRIT(s) = mean(err_esprit(~isnan(err_esprit)));
end

% === Plot RMSE vs. SNR ===
figure;
plot(snrVals, RMSE_MUSIC, '-ob', 'LineWidth', 2, 'DisplayName', 'MUSIC'); hold on;
plot(snrVals, RMSE_ESPRIT, '-sr', 'LineWidth', 2, 'DisplayName', 'ESPRIT');

% Highlight invalid ESPRIT points
nanIdx = isnan(RMSE_ESPRIT);
if any(nanIdx)
    scatter(snrVals(nanIdx), zeros(1,sum(nanIdx)), 100, 'rx', 'LineWidth', 2, ...
        'DisplayName', 'ESPRIT failed');
end

xlabel('SNR (dB)');
ylabel('RMSE (degrees)');
title('DOA Estimation RMSE vs. SNR');
legend('show', 'Location', 'northeast');
grid on;
