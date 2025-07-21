function [estAngle, x_power, MSEangle_dB] = SBLMC(R,Psi,Xi,zeta,K,M,N,P,U,isUpdateMC,isUpdateNu,isShowFig,targetAngle,Gamma)

%   SBLMC algorithm for DOA estimation in MIMO radar with unknown mutual coupling
%
%	@input param
%   R: received signal
%   Psi: dictionary matrix
%   Xi: first order of dictionary matrix
%   zeta: discretized angles
%   K: number of targets
%   M: number of transmitting antennas
%   N: number of receiving antennas
%   P: number of pulses
%   U: number of discretized angles
%   isUpdateMC: is update the mutual coupling vectors in processes
%   isUpdateNu: is update the off-grid angle vector in processes
%   isShowFig: show the figure in processes
%   targetAngle: target angles as ref DOAs
%   Gamma: sparse matrix as ref 
%
%	@output param
%   estAngle: estimated target DOAs
%   x_power: the spatial spectrum
%   MSEangle_dB: the DOA estimation error in each iterative
%
%	Reference: Off-Grid DOA Estimation Using Sparse Bayesian Learning in MIMO Radar With Unknown Mutual Coupling
%				Author: Peng Chen (chenpengseu@seu.edu.cn)
%				Date: 2018-04-12

vec = @(MAT) MAT(:);

%% SBLMC algorithm
% initial 
alphan_hat = mean(var(R));
cT_hat = [1; zeros(M-1,1)];
cR_hat = [1; zeros(N-1,1)];

nu_hat = zeros(U, 1);
beta_hat = 1./mean(abs(Psi(:,1:M*N:end)'*R), 2);

% beta_hat = mean(abs(Psi'*R), 2);
varthetaT_hat = [1;zeros(M-1,1)];
varthetaR_hat = [1;zeros(N-1,1)];

% a = 1e-2; 
a = 1e-4; 
b = a;
c = a; d =a; 
e1 = a; f1 = a;
e2 = a; f2= a;

tol = 1e-3;
MSEangle_dB = [];
for iter = 1 : 1 : 1e3
    iter
    % update mu_p Sigma_x
    c_hat = kron(cR_hat, cT_hat);
    Upsilon_hat = Psi+Xi*kron(diag(nu_hat), eye(M*N));
    % Uc_hat = Upsilon_hat * kron(c_hat, eye(U));
    Uc_hat = Upsilon_hat * kron(eye(U), c_hat);
    SigmaX_hat = inv(alphan_hat * Uc_hat' * Uc_hat + diag(beta_hat));
    MU_hat = alphan_hat * SigmaX_hat * Uc_hat' * R;

    x_power = real(diag(SigmaX_hat)) + mean(abs(MU_hat).^2, 2);

    if isShowFig
        if exist('FH1')
            delete(FH1);
            delete(FH2);
        end
        [sortAngle, sortIndex] = sort(zeta+nu_hat, 'ascend');
        figure(1), FH1 = semilogy(rad2deg(sortAngle), x_power(sortIndex), 'b', 'LineWidth', 2); 
        figure(1), hold on; 
        figure(1), grid on;
        figure(1), FH2 = stem(rad2deg(targetAngle), max(x_power) * ones(K,1), '*', 'r', 'LineWidth', 2, 'MarkerSize', 10);
        legend('Estimated spatial spectrum', 'Target angles');
        drawnow;
    end


    % update beta
    beta_hat_last = beta_hat;
    % beta_hat = (P+c-1)./(d + P*real(diag(SigmaX_hat))+ sum(abs(MU_hat).^2,2));
    beta_hat = (P - P*beta_hat.*real(diag(SigmaX_hat)))./ (d+sum(abs(MU_hat).^2,2));

    if isUpdateMC
        % update ct
        for idxk1 = 1 : 1 : U
            for idxk2 = 1 : 1 : U
                Upsilon_hat1=Upsilon_hat(:, M*N*(idxk1-1)+1:M*N*idxk1);
                Upsilon_hat2=Upsilon_hat(:, M*N*(idxk2-1)+1:M*N*idxk2);
                UpSumTemp = Upsilon_hat1'*Upsilon_hat2 * SigmaX_hat(idxk2, idxk1);
                if idxk1 == 1 && idxk2 == 1
                    UpSum = UpSumTemp;
                else
                    UpSum = UpSum + UpSumTemp;
                end
            end
        end
        for idxP = 1 : 1 : P
            GT = [];
            for idxM = 1 : 1 : M
                eM = zeros(M, 1);
                eM(idxM) = 1;
                gTtemp = kron(cR_hat, eM);
                if idxM == 1
                    GT = zeros(length(gTtemp), M);
                end
                GT(:, idxM) = gTtemp;
                % GT = [GT, gTtemp];
            end
            TT = [];
            for idxM = 1 : 1 : M
                eM = zeros(M, 1);
                eM(idxM) = 1;
                tTtemp = kron(kron(MU_hat(:, idxP), cR_hat), eM);
                if idxM == 1
                    TT = zeros(length(tTtemp), M);
                end
                TT(:, idxM) = tTtemp;
                % TT = [TT, tTtemp];
            end

            zTtemp = alphan_hat * TT' * Upsilon_hat'*R(:, idxP);
            if idxP == 1
                zT = zTtemp;
            else
                zT = zT + zTtemp;
            end

            UpTTsumTemp = alphan_hat * TT' * Upsilon_hat' * Upsilon_hat * kron(kron(MU_hat(:, idxP), cR_hat), eye(M));
            if idxP == 1
                UpTTsum = UpTTsumTemp;
            else
                UpTTsum = UpTTsum + UpTTsumTemp;
            end
        end
        HT = alphan_hat * P * GT' * UpSum' * kron(cR_hat, eye(M)) + diag(varthetaT_hat) + UpTTsum;
        cT_hat = vec(inv(HT) * zT);

        % update varthetaT
        % varthetaT_hat = (1+e1) ./ (f1 + abs(cT_hat).^2);
        varthetaT_hat = 1 ./ (f1 + abs(cT_hat).^2);

        % update cr
        for idxk1 = 1 : 1 : U
            for idxk2 = 1 : 1 : U
                Upsilon_hat1=Upsilon_hat(:, M*N*(idxk1-1)+1:M*N*idxk1);
                Upsilon_hat2=Upsilon_hat(:, M*N*(idxk2-1)+1:M*N*idxk2);
                UpSumTemp = Upsilon_hat1'*Upsilon_hat2 * SigmaX_hat(idxk2, idxk1);
                if idxk1 == 1 && idxk2 == 1
                    UpSum = UpSumTemp;
                else
                    UpSum = UpSum + UpSumTemp;
                end
            end
        end
        for idxP = 1 : 1 : P
            GR = [];
            for idxN = 1 : 1 : N
                eN = zeros(N, 1);
                eN(idxN) = 1;
                gRtemp = kron(eN, cT_hat);
                if idxN == 1
                    GR = zeros(length(gRtemp), N);
                end
                GR(:, idxN) = gRtemp;
                % GR = [GR, gRtemp];
            end
            TR = [];
            for idxN = 1 : 1 : N
                eN = zeros(N, 1);
                eN(idxN) = 1;
                tRtemp = kron(kron(MU_hat(:, idxP), eN), cT_hat);
                % TR = [TR, tRtemp];
                if idxN == 1
                    TR = zeros(length(tRtemp), N);
                end
                TR(:, idxN) = tRtemp;
            end

            zRtemp = alphan_hat * TR' * Upsilon_hat'*R(:, idxP);
            if idxP == 1
                zR = zRtemp;
            else
                zR = zR + zRtemp;
            end

            UpTRsumTemp = alphan_hat * TR' * Upsilon_hat' * Upsilon_hat * kron(kron(MU_hat(:, idxP), eye(N)), cT_hat);
            if idxP == 1
                UpTRsum = UpTRsumTemp;
            else
                UpTRsum = UpTRsum + UpTRsumTemp;
            end
        end
        HR = alphan_hat * P * GR' * UpSum' * kron(eye(N), cT_hat) + diag(varthetaR_hat) + UpTRsum;
        cR_hat = vec(inv(HR) * zR);

        % update varthetaR
        % varthetaR_hat = (1+e2) ./ (f2 + abs(cR_hat).^2);
        varthetaR_hat = 1 ./ (f2 + abs(cR_hat).^2);
    end

    % update nu
    if isUpdateNu
        H = zeros(U, U);
        for idxU = 1 : 1 : U
            for idxM = 1 : 1 : U
                XiM = Xi(:, M*N*(idxM-1)+1:M*N*idxM);
                XiU = Xi(:, M*N*(idxU-1)+1:M*N*idxU);
                sumMu = vec(MU_hat(idxM, :))' * vec(MU_hat(idxU, :));
                H(idxU, idxM) = real((P*SigmaX_hat(idxU, idxM) + sumMu) * c_hat'*XiM'*XiU*c_hat);
            end
        end
        z = zeros(U, 1);
        for idxU = 1 : 1 : U
            XiU = Xi(:, M*N*(idxU-1)+1:M*N*idxU);
            sumP = 0;
            for idxP = 1 : 1 : P
                sumP = sumP + real((R(:, idxP) - Psi*kron(MU_hat(:, idxP), c_hat))' * XiU*MU_hat(idxU,idxP)*c_hat);
            end
            sumM = 0;
            for idxM = 1 : 1 : U
                PsiM = Psi(:, M*N*(idxM-1)+1:M*N*idxM);
                sumM = sumM + real(P * SigmaX_hat(idxU, idxM) * c_hat'*PsiM'*XiU*c_hat);
            end
            z(idxU) = sumP - sumM;
        end
        nu_hat = inv(H) * z;
    end

    % update alphan
    UICR = Upsilon_hat * kron(eye(U), c_hat);
    % alphan_hat = (M*N*P+ a - 1) / (b+norm(R-Upsilon_hat*kron(MU_hat, c_hat),'fro')^2+P*real(trace(UICR' * UICR * SigmaX_hat)));
    alphan_hat = (M*N*P - 1-P*real(trace(UICR' * UICR * SigmaX_hat)) * alphan_hat) / (b+norm(R-Upsilon_hat*kron(MU_hat, c_hat),'fro')^2);



    if norm(beta_hat - beta_hat_last)/norm(beta_hat_last) < tol
        break;
    end

    % peak position
    peakIndex = [];
    for idxPower = 2 : 1 : length(x_power) - 1
        if x_power(idxPower) >= x_power(idxPower-1) && x_power(idxPower) >= x_power(idxPower+1)
            peakIndex = [peakIndex; idxPower];
        end
    end
    if length(peakIndex) < K
        [~, sortIndex] = sort(x_power, 'descend');
        for idxSort = 1 : 1 : length(sortIndex)
            if ~ismember(sortIndex(idxSort), peakIndex)
                peakIndex = [peakIndex; idxPower];
            end
            if length(peakIndex) >= K
                break;
            end
        end
    end

    % target angle
    [~, sortIndex] = sort(x_power(peakIndex), 'descend');
    angleIndexEst = peakIndex(sortIndex(1:K,1));
    estAngle = zeta(angleIndexEst)+nu_hat(angleIndexEst);
    allPerm = perms(1:K);
    [~,minIndex] = min(sum(abs(bsxfun(@minus, estAngle(allPerm).', targetAngle)).^2));
    estAngle = estAngle(allPerm(minIndex,:));

    % MSE angle
    MSEangle_dB = [MSEangle_dB; 10*log10(norm(estAngle - targetAngle)^2)];
end

