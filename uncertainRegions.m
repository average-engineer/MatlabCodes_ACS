%% Settings
clearvars
close all
clc
format short

%% Uncertain Parameter Ranges
% Controller Gain
kmin = 2;
kmax = 3;
k = [kmin:0.1:kmax];
% Time Delay
thetaMin = 2;
thetaMax = 3;
theta = [thetaMin:0.1:thetaMax];
% Time Constant
tauMin = 2;
tauMax = 3;
tau = [tauMin:0.1:tauMax];
%% Plant Family
s = tf('s');
count = 1;
for ii = 1:length(k)
    for jj = 1:length(theta)
        for kk = 1:length(tau)
            Gp(count) = (k(ii)/(tau(kk)*s + 1))*exp(-theta(jj)*s);
            count = count + 1;
        end
    end
end

%% Plotting the actual uncertainty region
point = 0.2*1i; % selected point
res = evalfr(Gp,point);

figure
plot(res,'*')
grid on
xlabel('Re')
ylabel('Im')