%% Settings
clearvars
close all
clc
format short
s = tf('s');
%% Nominal Plant
Gnom = 3*(1-2*s)/((5*s+1)*(10*s+1));
%% Controller
Kc1 = 1.13; % Ziegler Design
Kc2 = 0.31;
Kc = [Kc1;Kc2];
pI = (12.7*s+1)/(12.7*s); % PI Controller
%% Closed Loop System
L1 = Gnom*Kc1*pI;
L2 = Gnom*Kc2*pI;
T1 = L1/(1+L1);
T2 = L2/(1+L2);
%% Multiplicative Uncertainty Weight
wM = (10*s + 0.33)/(((10*s)/5.25) + 1);

%% Frequency response
freq = logspace(-3,2,100);
for kk = 1:length(freq)
    wMFR(kk) = abs(evalfr(1/wM,freq(kk)*1i)); % Inverse Weight Response
    T1FR(kk) = abs(evalfr(T1,freq(kk)*1i));
    T2FR(kk) = abs(evalfr(T2,freq(kk)*1i));
end
figure
semilogx(freq,20*log10(wMFR),'--','color','k','DisplayName','Upper Bound on T for RS')
hold on
semilogx(freq,20*log10(T1FR),'DisplayName','T1')
semilogx(freq,20*log10(T2FR),'DisplayName','T2')
grid on
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
legend
