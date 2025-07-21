clear all; format long; warning off;
% random seed
rng(222); 

vec = @(MAT) MAT(:);
% trans antennas
M = 10;
% recv antennas
N = 5;
% pulse num
P = 100;
% target num
% K = 3;
K = 3;

% antenna space 
dT = 0.5;
dR = dT;

% steering vectors
steerT = @(angleT) exp(1j*2*pi*dT* [0 : 1 : M-1].' * sin(angleT(:).'));
steerR = @(angleR) exp(1j*2*pi*dR* [0 : 1 : N-1].' * sin(angleR(:).'));

% grids
delta = deg2rad(2);
zeta = deg2rad([-80:rad2deg(delta):80].');
U = length(zeta);

% true DOA
angleIndex = vec(randperm(U-1, K));
targetAngleGrid = vec(zeta(angleIndex));

isOffGrid = true;
minSpace = deg2rad(20);
while true
	if isOffGrid
		targetAngleDelta = rand(K,1) * delta - delta / 2;
	else
		targetAngleDelta = zeros(K,1);
	end
	targetAngle = targetAngleGrid + targetAngleDelta;
	angSort = sort(targetAngle, 'ascend');
	if min(abs(angSort(2:end) - angSort(1:end-1))) >= minSpace
		break;
	end
end
theta = targetAngle;


% true target scattering coefficients
Gamma = (ones(K, P)+0.1*unifrnd(-0.5,0.5,K,P)) .* exp(1j*2*pi*unifrnd(0,1,K,P));

% mutual coupling matrix
mutualRef_dB = -5;
mutualVecT = 10.^((mutualRef_dB * [1 : 1 : M - 1].') / 20);
mutualVecT = mutualVecT .* (1 + unifrnd(-0.05,0.05,M-1, 1));
mutualVecT = sort(mutualVecT, 'descend');
mutualVecT = mutualVecT .* exp(1j * 2 * pi * rand(M-1, 1));
cT = [1; mutualVecT];
mutualVecR = 10.^((mutualRef_dB * [1 : 1 : N - 1].') / 20);
mutualVecR = mutualVecR .* (1 + unifrnd(-0.05, 0.05, N-1, 1));
mutualVecR = sort(mutualVecR, 'descend');
mutualVecR = mutualVecR .* exp(1j * 2 * pi * rand(N-1, 1));
cR = [1; mutualVecR];

% generate mutual coupling matrix
Ct = zeros(M, M);
for idxM1 = 1 : 1 : M
    for idxM2 = 1 : 1 : M
        Ct(idxM1, idxM2) = cT(abs(idxM1-idxM2) + 1);
    end
end
Cr = zeros(N, N);
for idxN1 = 1 : 1 : N
    for idxN2 = 1 : 1 : N
        Cr(idxN1, idxN2) = cR(abs(idxN1-idxN2) + 1);
    end
end

% echo signal
for idxTarget = 1 : 1 : K
    Dtemp = kron(steerR(theta(idxTarget)), steerT(theta(idxTarget)));
    if idxTarget == 1
        D = zeros(length(Dtemp), K);
    end
    D(:, idxTarget) = Dtemp;
end
echoSignal = kron(Cr, Ct) * D * Gamma;

% add awgn
SNR_dB = 20;
noiseVariance = echoSignal(:)' * echoSignal(:) / length(echoSignal(:)) / db2pow(SNR_dB);
noise =  sqrt(noiseVariance / 2) * (randn(size(echoSignal)) + 1j * randn(size(echoSignal)));
R = echoSignal + noise;


% dictionary matrix
Psi = [];
Xi = [];
for idxU = 1 : 1 : U
    % Qa matrix
    aTemp = steerT(zeta(idxU));
    Qa1 = zeros(M);
    Qa2 = zeros(M);
    partialQa1 = zeros(M);
    partialQa2 = zeros(M);
    for p = 1 : 1 : M
        for q = 1 : 1 : M
            Qa1(p,q) = 0;
            if p+q<=M+1
                Qa1(p,q) = aTemp(p+q-1);
            end
            Qa2(p,q) = 0;
            if p>=q && q>=2
                Qa2(p,q) = aTemp(p-q+1);
            end

            partialQa1(p,q) = 0;
            if p+q<=M+1
                partialQa1(p,q) = 1j*2*pi*(p+q-2)*dT*cos(zeta(idxU))*aTemp(p+q-1);
            end
            partialQa2(p,q) = 0;
            if p>=q && q>=2
                partialQa2(p,q) = 1j*2*pi*(p-q)*dT*cos(zeta(idxU))*aTemp(p-q+1);
            end

        end
    end
    Qa = Qa1+Qa2;
    partialQa = partialQa1+partialQa2;

    % Qb matrix
    bTemp = steerR(zeta(idxU));
    Qb1 = zeros(N);
    Qb2 = zeros(N);
    partialQb1 = zeros(N);
    partialQb2 = zeros(N);
    for p = 1 : 1 : N
        for q = 1 : 1 : N
            Qb1(p,q) = 0;
            if p+q<=N+1
                Qb1(p,q) = bTemp(p+q-1);
            end
            Qb2(p,q) = 0;
            if p>=q && q>=2
                Qb2(p,q) = bTemp(p-q+1);
            end

            partialQb1(p,q) = 0;
            if p+q<=N+1
                partialQb1(p,q) = 1j*2*pi*(p+q-2)*dR*cos(zeta(idxU))*bTemp(p+q-1);
            end
            partialQb2(p,q) = 0;
            if p>=q && q>=2
                partialQb2(p,q) = 1j*2*pi*(p-q)*dR*cos(zeta(idxU))*bTemp(p-q+1);
            end
        end
    end
    Qb = Qb1+Qb2;
    partialQb = partialQb1+partialQb2;

    % Phi matrix
    Phi = kron(Qb, Qa);
    partialPhi = kron(Qb, partialQa) + kron(partialQb, Qa);
    Psi = [Psi, Phi];
    Xi = [Xi, partialPhi];
end
% mutual coherence

% echo signal
for idx = 1 : 1 : length(zeta)
    Dtemp = kron(steerR(zeta(idx)), steerT(zeta(idx)));
    % Dtemp = steerR(zeta(idx));
    if idx == 1
        DD = zeros(length(Dtemp), length(zeta));
    end
    DD(:, idx) = Dtemp;
end




isUpdateMC=true;
isUpdateNu=true;
isShowFig=true;
% SBLMC algorithm
[estAngle, x_power, MSEangle_dB] = SBLMC(R,Psi,Xi,zeta,K,M,N,P,U,isUpdateMC, isUpdateNu,isShowFig,targetAngle,Gamma);

% MUSIC algorithm
% [estAngle, MUSICpower, MSEangleMUSIC_dB] = MUSIC(R,K,steerT,steerR,isShowFig,targetAngle);

