%% Settings
clearvars
close all
clc
format short
s = tf('s');
%% Nominal Loop Gain
L1 = 0.5/s;
L2 = (0.5/s)*((1-s)/(1+s));
%% Closed Loop Sensitivity TF
S1 = 1/(1+L1);
S2 = 1/(1+L2);
%% Performance Weight 
wP = 0.25 + (0.1/s);
%% Uncertainty Weight
wU = 0.85*s/(s+1);
%% Frequency Repsonse
freq = logspace(-3,5,100);
for kk = 1:length(freq)
    S1FR(kk) = abs(evalfr(S1,freq(kk)));
    S2FR(kk) = abs(evalfr(S2,freq(kk)));
    wPFR(kk) = abs(evalfr(wP,freq(kk)));
    wUFR(kk) = abs(evalfr(wU,freq(kk)));
end
S1peak = max(S1FR)
S2peak = max(S2FR)
%% Performance
figure
loglog(freq,1./wPFR,'--','color','k','DisplayName','Nominal Performance Bound')
hold on
loglog(freq,1./((wPFR + wUFR)),'--','color','r','DisplayName','Robust Performance Bound')
loglog(freq,S1FR,'DisplayName','Nominal Loop Gain 1')
loglog(freq,S2FR,'DisplayName','Nominal Loop Gain 2')
grid on
xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('Performance')
legend
