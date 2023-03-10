%% Author: Ashutosh Mukherjee
% Co-Author: Raj Khamkar
% Last Date of modification: 26.02.2023
%% Settings
clearvars
close all
clc
format long
addpath("matlabtikz")
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
%% System Reduction
[hsv_stab, hsv_unstab] = hankelsv(sys, 'ncf', 'log'); % Small value correspond to states that do not contribute much to future state trajectory
sys_red = reduce(sys, 4, 'errortype', 'ncf'); % Reduce model to 4th-order
Gred = minreal(tf(sys_red));
%% Checking pre-requisits for Doubly Co-prime Factorization
n = 4; % # states
m = size(Gred,2); % # inputs
p = size(Gred,1); % # outputs
if rank(ctrb(sys_red)) == n && rank(obsv(sys_red)) == n
    fprintf("Pre-requisites met\n");
else
    fprintf("Pre-requisites not met\n");
end

%% F matrix
% Af = A + B*F needs to be stable
desEigs = [-1,-1,-2,-2]'; % Desired eigenvalues of Af
F = place(sys_red.A,-sys_red.B,desEigs);
Af = sys_red.A + sys_red.B*F;
Cf = sys_red.C + sys_red.D*F;
%% H matrix
% Ah = A + H*C needs to be stable
H = place(sys_red.A',-sys_red.C',desEigs)';
Ah = sys_red.A + H*sys_red.C;
Bh = sys_red.B + H*sys_red.D;
%% Doubly Co-Prime Factorization
N = tf(ss(Af,sys_red.B,Cf,sys_red.D));
M = tf(ss(Af,sys_red.B,F,eye(m,m)));
X = tf(ss(Af,-H,Cf,eye(m,m)));
Y = tf(ss(Af,-H,F,zeros(m,m)));

Ntilda = tf(ss(Ah,Bh,sys_red.C,sys_red.D));
Mtilda = tf(ss(Ah,H,sys_red.C,eye(m,m)));
Xtilda = tf(ss(Ah,-Bh,F,eye(m,m)));
Ytilda = tf(ss(Ah,-H,F,zeros(m,m)));
%% Verifying
ver = [Xtilda,-Ytilda;-Ntilda,Mtilda]*[M,Y;N,X];
freq = logspace(-2,3,100);
for kk = 1:length(freq)
    verFR(:,:,kk) = evalfr(ver,freq(kk)*1i); % to ensure that identity matrix is obtained
end
%% Youla Controller 
% Q = [tf(1,[1,2]),-tf(2,[1,2]);-tf(2,[1,2]),tf(2,[1,3])]; % Stable, real-rational proper TF
Q = [tf(1,[1,1]),0;0,tf(1,[1,1])];
% Q = tf(1,[1 8 15])*eye(size(Gred));
K = (Y - M*Q)*inv(X - N*Q);
% K = inv(Xtilda-Q*Ntilda)*(Ytilda-Q*Mtilda);
K = minreal(K);
%% Closed Loop Transfer Functions
G = N*inv(M);
L = G*K;
% S = inv(eye(size(L)) - L);
S = Mtilda*(X-Ntilda*Q);
T = Ntilda*(Mtilda*Q-Y);
%% Closed Loop Time Response
% Arbritary Input
t = linspace(0,10,100);
u1 = 2.*exp(-t);
u2 = zeros(size(t));
u = [u1;u2];
figure
subplot(2,1,1)
y1 = lsim(S,u,t);
plot(t,y1(:,1),'color','k','linewidth',2,'DisplayName','S y1')
hold on
plot(t,y1(:,2),'color','r','linewidth',2,'DisplayName','S y2')
grid on
legend
xlabel('Time (s)')
ylabel('Amplitude (-)')
subplot(2,1,2)
y2 = lsim(T,u,t);
plot(t,y2(:,1),'color','k','linewidth',2,'DisplayName','T y1')
hold on
plot(t,y2(:,2),'color','r','linewidth',2,'DisplayName','T y2')
grid on
legend
xlabel('Time (s)')
ylabel('Amplitude (-)')
grid on
legend
matlab2tikz();
%% Singular Value Plots
for kk = 1:length(freq)
    [~,eeS,~] = svd(evalfr(S,freq(kk)*1i));
    [~,eeT,~] = svd(evalfr(T,freq(kk)*1i));
    % Upper Singular Values
    usvS(kk) = eeS(1,1);
    usvT(kk) = eeT(1,1);
    % Lower Singular Values
    lsvS(kk) = eeS(end,end);
    lsvT(kk) = eeT(end,end);
end

figure
subplot(2,1,1)
semilogx(freq,usvS,'--','color','k','linewidth',2,'DisplayName','Senstivity TF')
hold on
semilogx(freq,usvT,'linewidth',2,'DisplayName','Complimentary Sensitivity TF')
grid on
xlabel('Frequency (rad/s)')
ylabel('Upper Singular Value (-)')
legend
subplot(2,1,2)
semilogx(freq,lsvS,'--','color','k','linewidth',2,'DisplayName','Senstivity TF')
hold on
semilogx(freq,lsvT,'linewidth',2,'DisplayName','Complimentary Sensitivity TF')
grid on
xlabel('Frequency (rad/s)')
ylabel('Lower Singular Value (-)')
legend
% matlab2tikz()