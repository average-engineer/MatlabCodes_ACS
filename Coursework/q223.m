%% Author: Ashutosh Mukherjee
% Co-Author: Raj Khamkar
% Last Date of modification: 26.02.2023
%% Settings
clearvars
close all
clc
format short
addpath("matlabtikz")
%% Plant
G = tf([1 -1],[1 1 -4 -4]);
%% Verification of Bezout Identity
N = tf(zpk(1,[-1 -1 -1],1));
M = tf(zpk([-2 2],[-1 -1],1));
U = tf([61 121],[3 3]);
V = tf([3 12 -31],[3 6 3]);
bzIden = U*N + V*M
%% Youla Stabilizing Controller
Q = tf(1,[1 3]);
K = (U + Q*M)/(V - Q*N);
%% Closed Loop TFs
Tf1 = minreal(inv(eye(size(G)) + (K*G)));
Tf2 = minreal(K*inv(eye(size(G)) + (G*K)));
Tf3 = minreal(G*Tf1);
Tf4 = minreal(inv(eye(size(G)) + (G*K)));
%% Senstivity and Complimentary Sensitivity Transfer Function
S = M*(V-Q*N);
T = N*(U+Q*M);
%% Closed loop response
t = linspace(0,20,100);
u = 2.*exp(-t);
figure
subplot(2,2,1)
x = lsim(Tf1,u,t);
plot(t,x,'linewidth',2,'DisplayName','TF1')
xlabel('Time(s)')
ylabel('Amplitude(-)')
legend
grid on
subplot(2,2,2)
x = lsim(Tf2,u,t);
plot(t,x,'linewidth',2,'DisplayName','TF2')
xlabel('Time(s)')
ylabel('Amplitude(-)')
legend
grid on
subplot(2,2,3)
x = lsim(Tf3,u,t);
plot(t,x,'linewidth',2,'DisplayName','TF3')
xlabel('Time(s)')
ylabel('Amplitude(-)')
legend
grid on
subplot(2,2,4)
x = lsim(Tf4,u,t);
plot(t,x,'linewidth',2,'DisplayName','TF4')
xlabel('Time(s)')
ylabel('Amplitude(-)')
legend
grid on
matlab2tikz();

figure
xs = lsim(S,u,t);
plot(t,xs,'Color','r','linewidth',2,'DisplayName','Sensitivity Transfer Function')
hold on
xt = lsim(T,u,t);
plot(t,xt,'Color','k','linewidth',2,'DisplayName','Complimentary Sensitivity Transfer Function')
grid on
xlabel('Time(s)')
ylabel('Amplitude(-)')
legend
%matlab2tikz();