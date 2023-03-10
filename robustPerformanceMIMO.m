%% Settings
clearvars
close all
clc
format short
s = tf('s');
freq = logspace(-3,2,100);
%% Nominal Plant
% G = (1/(75*s + 1))*[87.8,-86.4;108.2,-109.6]; % Distillation Process
G = (1/(75*s + 1))*[87.8,1.4;-108.2,-1.4]; % Distillation Process

% Condition Number
for kk = 1:length(freq)
    cN(kk) = cond(evalfr(G,freq(kk)));
end
%% Controller
K = (0.7/s)*inv(G); % Inverse Based Controller
%% Weights
wP = (0.5*s + 0.05)/(s); % Performance Weight
wI = (s + 0.2)/(0.5*s + 1); % Input multiplicative uncertainty weight
%% Nominal Stability
L = G*K; % Nominal Loop Gain
S = inv(eye(size(L)) + L); % Output Sensitivity TF
TI = K*G*inv(eye(size(L)) + (K*G)); % Input Complimentary Sensitivity TF
% Step Responses
% figure
% subplot(2,1,1)
% step(TI)
% title('TI')
% grid on
% subplot(2,1,2)
% step(S)
% title('S')
% grid on
%% Nominal Performance
for kk = 1:length(freq)
    nom = evalfr(wP*S,freq(kk)*1i);
    % Upper singular value of S at each frequency
    [~,ee,~] = svd(nom);
    nomFR(kk) = ee(1,1);
end

%% Robust Stability
% Structured Perturbation
% Input Perturbation = diag(del1,del2), del1,del2 : complex scalars
BlockStructure1 = [1,0;1,0];

% Unstructured Perturbation
BlockStructure2 = [2,2];
for kk = 1:length(freq)
    M = wI*TI;
    MFR = evalfr(M,freq(kk)*1i);
    bounds1 = mussv(MFR,BlockStructure1);
    bounds2 = mussv(MFR,BlockStructure2);
    SSVRS1(kk) = max(bounds1); % Upper Bound selected as it is convex optimization -> More reliable
    SSVRS2(kk) = max(bounds2); 
end

%% Robust Performance
% Perturbation Structure 1
% Perturbation Block 1: del1 (Complex Scalar)
% Perturbation Block 2: del2 (Complex Scalar)
% Perturbation Block 3: Peformance Full Complex Matrix Perturbation of 2x2
BlockStructure1 = [1,0;1,0;2,2];

% Perturbation Structure 2
% Perturbation Block 1: Full Complex Matrix Perturbation of 2x2
% Perturbation Block 2: Peformance Full Complex Matrix Perturbation of 2x2
BlockStructure2 = [2,2;2,2];

% Closing the control loop using lower LFT
N = [wI*TI,wI*K*S;wP*S*G,wP*S];
N = minreal(N);

for kk = 1:length(freq)
    NFR = evalfr(N,freq(kk)*1i);
    bounds1 = mussv(NFR,BlockStructure1);
    bounds2 = mussv(NFR,BlockStructure2);
    SSVRP1(kk) = max(bounds1); % Upper Bound selected as it is convex optimization -> More reliable
    SSVRP2(kk) = max(bounds2);
end

%% Frequency response
figure
semilogx(freq,nomFR,'color','b','DisplayName','NP')
hold on
semilogx(freq,SSVRS1,'color','r','DisplayName','RS (Structured Perturbation)')
semilogx(freq,SSVRS2,'--','color','r','DisplayName','RS (Unstructured Perturbation)')
semilogx(freq,SSVRP1,'color','k','DisplayName','RP (Structured Perturbation)')
semilogx(freq,SSVRP2,'--','color','k','DisplayName','RP (Unstructured Perturbation)')
grid on
xlabel('Frequeny (rad/s)')
legend