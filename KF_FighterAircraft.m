%% Settings
clearvars
close all
clc
format short
%% Linearised Dynamics
A = [-1.7,50,260;
    0.22,-1.4,-32;
    0,0,-12];
B = [-272;0;14];
C = [1,0,0;
    0,1,0];
D = zeros(2,1);
F = [0.02,0.1;
    -0.0035,0.004;
    0,0]; % Process Noise Matrix
%% White Noise
V = F'*F; % Spectral density matrix of process noise
Z = 10*(C*C'); % Spectral density matrix of measurement noise
%% Kalman Filter
[Ke,Pe,eigKF] = lqe(A,F,C,V,Z)
%% Simulating white noise
randn('seed',0);
t = 0:0.01:10;
v = randn(size(t,2),2); % Process Noise
z = randn(size(t,2),2); % Measurement Noise
w = F*v' - Ke*z'; % Overall white noise in the system

figure
subplot(2,1,1)
plot(t,z(:,1))
grid on
xlabel('Time (s)')
title('Simulated Measurement White Noise')

subplot(2,1,2)
plot(t,v(:,1))
grid on
xlabel('Time (s)')
title('Simulated Process White Noise')
%% Estimation Error Dynamics
Ae = A - Ke*C;
Be = eye(3,3); % Input -> Overall White Noise in the system
Ce = eye(3,3); % Output is the estimation error vector itself
De = zeros(3,3);
sysE = ss(Ae,Be,Ce,De);
[e,T] = lsim(sysE,w',t);

figure
plot(T,e(:,1),'Linewidth',2,'DisplayName','e_1(t)')
hold on
plot(T,e(:,2),'Linewidth',2,'DisplayName','e_2(t)')
plot(T,e(:,3),'Linewidth',2,'DisplayName','e_3(t)')
grid on
xlabel('Time (s)')
title('Estimation Error Dynamics')
legend
%% Covariance of the computed estimation error
Pe1 = cov(e)