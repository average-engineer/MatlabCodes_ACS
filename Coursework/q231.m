%% Author: Ashutosh Mukherjee
% Co-Author: Raj Khamkar
% Last Date of modification: 26.02.2023
%% Settings
clearvars
close all
clc
format short
addpath("matlabtikz")
freq = logspace(-3,3,100);
uncertainty = "MitTD"; % OhneTD oder MitTD
%% Parametric Uncertainties
del = linspace(-1,1,21);
% Mass
mnom = 10;
mmax = 18;
mmin = 5;
bm = 0.385*mnom;
dm = -0.23;
for ii = 1:length(del)
    m(ii) = (mnom + bm*del(ii))/(1 + dm*del(ii));
end
% Spring Stiffness
knom = 1;
kmax = 1.8;
kmin = 0.8;
bk = -0.28*knom;
dk = -0.6;
delk = linspace(-1,1,21);
for ii = 1:length(del)
    k(ii) = (knom + bk*del(ii))/(1 + dk*del(ii));
end
% Output Gain
hnom = 2.5;
bh = -0.28*hnom;
dh = -0.6;
for ii = 1:length(del)
    h(ii) = (hnom + bh*del(ii))/(1 + dh*del(ii));
end
% Damping
cnom = 10;
cmax = 14;
cmin = 6;
rc = (cmax-cmin)/(cmax + cmin);
for ii = 1:length(del)
    c(ii) = cnom*(1 + rc*del(ii));
end

mp = ureal('m',mnom,'Percentage',[-50,80]);
cp = ureal('c',cnom,'Percentage',[-40,40]);
kp = ureal('k',knom,'Percentage',[-20,80]);
hp = ureal('h',hnom,'Percentage',[-20,80]);
%% Dynamic Uncertainty
% Time Delay
taunom = 0.1;
taumax = 1.3*taunom;
% Uncertainty weight (Lundstörm)
wTau = tf([2*taumax 0],[taumax 2]);
% Complex Perturbation for time delay
deltau = ultidyn('delT',[1,1]);

%% Weights
W1 = tf([2 0],[1 0.1127]); % Output Weight Weight
W2 = tf([0.7 1],[0.07 1]); % Actuator Effort Weight
W2 = 0.3*W2;
%% PI Controller Design for Nominal System
K0 = 5.3547*tf([10 1],[10 0]); % Lag Compensator (designed so as to keep the gain cross-over frequency at 1 rad/s)
K = W2*K0; % W2 acting as a Lead Compensator
Lnom = minreal(W1*hnom*K*tf(1,[mnom cnom knom]));

%% Nominal Stability
clTF1 = W1/(1 + Lnom); % closed loop TF from disturbance to output
for kk = 1:length(freq)
    LnomFR(kk) = abs(evalfr(Lnom,freq(kk)*1i));
end

figure
semilogx(freq,LnomFR,'--','color','k','DisplayName','Lead-Lag Compensated Loop Gain')
xlabel('Frequency (rad/s)')
grid on
legend

figure
step(clTF1)
grid on

switch uncertainty
    case "OhneTD"
        n = 4; % number of uncertainties
        % Perturbation Block Structure
        % All 4 Perturbation Blocks: Real single perturbation
        BlockStructure = [-1,0;-1,0;-1,0;-1,0];
        % Manually obtained generalized plant
        P = [-bm/mnom,-1,-(bk - knom*dk),0,-knom,-cnom,0,1;
            zeros(1,5),rc*cnom,0,0;
            0,0,-dk,0,1,zeros(1,3);
            zeros(1,3),-dh,1,zeros(1,3);
            zeros(1,5),1,0,0;
            (dm/mnom) - (bm/(mnom^2)),-1/mnom,-(bk-knom*dk)/mnom,0,-knom/mnom,-cnom/mnom,0,1/mnom;
            zeros(1,3),W1*(bh-hnom*dh),W1*hnom,0,W1,0;
            zeros(1,3),-(bh-hnom*dh),-hnom,0,-1,0];
        % Uncertain Closed Loop
        Lp = W1*hp*K*tf(1,[mp cp kp]);
        clTFp = uss(W1/(1 + Lp));
    case "MitTD"
        n = 5; % number of uncertainties
        % Perturbation Block Structure
        % First 4 Perturbation Blocks: Real scalar perturbation
        % Last Perturbation Block: Complex scalar perturbation
        BlockStructure = [-1,0;-1,0;-1,0;-1,0;1,0];
        % Manually obtained generalized plant
        P = [-bm/mnom,-1,-(bk - knom*dk),0,1,-knom,-cnom,0,1;
            zeros(1,6),rc*cnom,0,0;
            0,0,-dk,0,0,1,zeros(1,3);
            zeros(1,3),-dh,0,1,zeros(1,3);zeros(1,8),wTau;
            zeros(1,6),1,0,0;
            (dm/mnom) - (bm/(mnom^2)),-1/mnom,-(bk-knom*dk)/mnom,0,1/mnom,-knom/mnom,-cnom/mnom,0,1/mnom;
            zeros(1,3),W1*(bh-hnom*dh),0,W1*hnom,0,W1,0;
            zeros(1,3),-(bh-hnom*dh),0,-hnom,0,-1,0];
        % Uncertain Closed Loop
        Lp = W1*hp*K*tf(1,[mp cp kp])*(1 + wTau*deltau);
        clTFp = uss(W1/(1 + Lp));
end


%% Closing the Lower Controller Loop
Nmanual = lft(P,K);
Mmanual = minreal(Nmanual(1:n,1:n));
%% Using Robust Control Toolbox
% Uncertain Closed Loop
[N,delta] = lftdata(clTFp); % taking out the uncertainty
N = minreal(tf(N));
N22 = minreal(N(end-2:end,end-2:end));
M = N(1:size(delta,1),1:size(delta,1));

%% Robust Stability
for kk = 1:length(freq)
    % Using upper singular value
    [~,ee,~] = svd(evalfr(Mmanual,freq(kk)*1i));
    [~,ee1,~] = svd(evalfr(M,freq(kk)*1i));
    MSVFR(kk) = ee(1,1);
    M1SVFR(kk) = ee1(1,1);
    % Using Structure Singular Value
    bounds = mussv(evalfr(Mmanual,freq(kk)*1i),BlockStructure);
    MSSVFR(kk) = max(bounds); % Using the upper bound (more reliable due to convex optimization)
    bounds1 = mussv(evalfr(M,freq(kk)*1i),BlockStructure);
    M1SSVFR(kk) = max(bounds1); % Using the upper bound (more reliable due to convex optimization)
end

figure
semilogx(freq,MSVFR,'--','color','k','linewidth',2,'DisplayName','Manually Computed M')
hold on
semilogx(freq,M1SVFR,'color','r','Linewidth',2,'DisplayName','M from Robust Control Toolbox')
xlabel('Frequency (rad/s)')
ylabel("Upper Singular Value")
grid on
legend


figure
semilogx(freq,M1SSVFR,'color','r','DisplayName','M from Robust Control Toolbox')
xlabel('Frequency (rad/s)')
ylabel("Structured Singular Value")
grid on
legend
matlab2tikz();

%% Nominal Performance
for kk = 1:length(freq)
    [~,ee,~] = svd(evalfr(N22,freq(kk)*1i));
    N22FR(kk) = ee(1,1);
end
figure
semilogx(freq,N22FR,'color','r','DisplayName','Upper Singular Value (N22)')
xlabel('Frequency (rad/s)')
ylabel("Upper Singular Value")
grid on
legend






