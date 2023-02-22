%% Settings
clearvars
close all
clc
format short
sympref('FloatingPointOutput',true); % format short equivalent for symolic math toolbox
syms s
%% Linearized Dynamics
ag = [-2.2567e-02  -3.6617e+01  -1.8897e+01  -3.2090e+01   3.2509e+00  -7.6257e-01
       9.2572e-05  -1.8997e+00   9.8312e-01  -7.2562e-04  -1.7080e-01  -4.9652e-03
       1.2338e-02   1.1720e+01  -2.6316e+00   8.7582e-04  -3.1604e+01   2.2396e+01
       0            0   1.0000e+00            0            0            0 
       0            0            0            0  -3.0000e+01            0 
       0            0            0            0            0  -3.0000e+01];
bg = [0     0 
      0     0 
      0     0 
      0     0 
     30     0 
      0    30];
cg = [0     1     0     0     0     0 
      0     0     0     1     0     0];
dg = [0     0 
      0     0];
sys = ss(ag, bg, cg, dg);

%% Reduce model order
[hsv_stab, hsv_unstab] = hankelsv(sys, 'ncf', 'log'); % Small value correspond to states that do not contribute much to future state trajectory
sys_red = reduce(sys, 4, 'errortype', 'ncf'); % Reduce model to 4th-order

% Rearrange to zpk-form
sys_zpk = zpk(sys_red);
[z, p, k1] = zpkdata(sys_zpk);

% Rounding up/down all zeros and poles and gains for simplicity
% z{1} = round(z{1}*100)/100;
% z{2} = round(z{2}*100)/100;
% z{3} = round(z{3}*100)/100;
% z{4} = round(z{4}*100)/100;
% k1   = round(k*100)/100;
% p{1} = round(p{1}*100)/100;
% p{2} = round(p{2}*100)/100;
% p{3} = round(p{3}*100)/100;
% p{4} = round(p{4}*100)/100;

%% Polynomial Matrix
% All poles from all G elements are same
N = [k1(1,1)*(s - z{1}(1))*(s - z{1}(2))*(s - z{1}(3)), k1(1,2)*(s - z{3}(1))*(s - z{3}(2))*(s - z{3}(3));
     k1(2,1)*(s - z{2}(1))*(s - z{2}(2))*(s - z{2}(3)), k1(2,2)*(s - z{4}(1))*(s - z{4}(2))*(s - z{4}(3))];
% Common Denominator Polynomial
d = (s - p{1}(1))*(s - p{1}(2))*(s - p{1}(3))*(s - p{1}(4));
 
%% Smith Form of Polynomial Matrix
S = smithForm(N);
%% Smith-McMillan Form of Reduced System
M = S/d;

%% Transfer Functions
Gfull = minreal(tf(sys));
Gred = minreal(tf(sys_red));
%% Singular Value Plots
freq = logspace(-3,3,100);
for kk = 1:length(freq)
    [~,eefull,~] = svd(evalfr(Gfull,freq(kk)*1i));
    [~,eered,~] = svd(evalfr(Gred,freq(kk)*1i));
    % Upper Singular Values
    usvFull(kk) = eefull(1,1);
    usvRed(kk) = eered(1,1);
    % Lower Singular Values
    lsvFull(kk) = eefull(end,end);
    lsvRed(kk) = eered(end,end);

    % Condition Number
    condFull(kk) = cond(evalfr(Gfull,freq(kk)*1i));
    condRed(kk) = cond(evalfr(Gred,freq(kk)*1i));
end

figure
subplot(2,1,1)
semilogx(freq,usvFull,'--','linewidth',2,'DisplayName','Full-Order Model')
hold on
semilogx(freq,usvRed,'linewidth',1,'DisplayName','Reduced-Order Model')
grid on
xlabel('Frequency (rad/s)')
title('Upper Singular Values')
legend
subplot(2,1,2)
semilogx(freq,lsvFull,'--','linewidth',2,'DisplayName','Full-Order Model')
hold on
semilogx(freq,lsvRed,'linewidth',1,'DisplayName','Reduced-Order Model')
grid on
xlabel('Frequency (rad/s)')
title('Lower Singular Values')
legend
figure
semilogx(freq,condFull,'--','linewidth',2,'DisplayName','Full-Order Model')
hold on
semilogx(freq,condRed,'linewidth',1,'DisplayName','Reduced-Order Model')
grid on
xlabel('Frequency (rad/s)')
title('Condition Number')
legend