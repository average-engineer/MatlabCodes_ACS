%% Settings
clearvars
close all
clc
format short
%% Linearized Longitudinal Dynamics of Flexible Bomber Aircraft
A = [0.4158,1.025,-0.00267,-0.0001106,-0.08021,0;
    -5.5,-0.8302,-0.06549,-0.0039,-5.115,0.809;
    0,0,0,1,0,0;
    -1040,-78.35,-34.83,-0.6214,-865.6,-631;
    0,0,0,0,-75,0;
    0,0,0,0,0,-100];
B = [zeros(4,2);75,0;0,100];
C = [-1491,-146.43,-40.2,-0.9412,-1285,-564.66;
    0,1,zeros(1,4)];
D = zeros(2,2);
F = B; % Process Noise Coefficient Matrix
n = size(A,1); % number of states
m = size(B,2); % number of inputs
l = size(C,1); % number of outputs
%% LQR Design
% Actuator Effort is 100 times more important than response decay
R = eye(m,m);
Q = 0.01*eye(n,n);
[K,P,lqrEig] = lqr(A,B,Q,R)
%% KF Design
V = 0.0007*(B'*B); % Process Noise Spectral Density Matrix
Z = C*C'; % Measurement Noise Spectral Density Matrix
[Ke,Pe,kfEig] = lqe(A,F,C,V,Z)
%% LQG
Alqg = A - B*K - Ke*C + Ke*D*K;
Blqg = Ke;
Clqg = -K;
Dlqg = zeros(2,2);
%% Closed-Loop System
sysP = ss(A,B,C,D); % Plant State Space
sysLQG = ss(Alqg,Blqg,Clqg,Dlqg); % LQG State Space
sysCL = feedback(sysP,sysLQG,+1); % Closed-Loop State Space
%% Closed-Loop Response
% CL State = [actual state;estimated state];
% Initial Condition
x0 = [0;0.1;zeros(10,1)];
[y,t,x] = initial(sysCL,x0); % x: [actual state;estimated state]
for ii = 1:length(t)
    u(ii,:) = -K*x(ii,7:12)'; % Actuator Effort = -K*(estimated state)
end

%% Full-State Feedback Response
x01 = [0;0.1;zeros(4,1)];
% Full State-Regulator (LQR) closed-loop system
sysLQR = ss(A-B*K,zeros(size(B)),C-D*K,zeros(size(D)));
[y1,t1,x1] = initial(sysLQR,x01); % x1: actual state vector
for ii = 1:length(t1)
    u1(ii,:) = -K*x1(ii,:)'; % Actuator Effort = -K*(estimated state)
end

%% LQG and LQR Comparison
figure
subplot(2,2,1)
plot(t,y(:,1),'DisplayName','LQG','LineWidth',2)
hold on
plot(t1,y1(:,1),'DisplayName','LQR','LineWidth',2)
legend
grid on
title("y_1(t)")

subplot(2,2,2)
plot(t,y(:,2),'DisplayName','LQG','LineWidth',2)
hold on
plot(t1,y1(:,2),'DisplayName','LQR','LineWidth',2)
legend
grid on
title("y_2(t)")

subplot(2,2,3)
plot(t,u(:,1),'DisplayName','LQG','LineWidth',2)
hold on
plot(t1,u1(:,1),'DisplayName','LQR','LineWidth',2)
legend
grid on
title("u_1(t)")

subplot(2,2,4)
plot(t,u(:,2),'DisplayName','LQG','LineWidth',2)
hold on
plot(t1,u1(:,2),'DisplayName','LQR','LineWidth',2)
legend
grid on
title("u_2(t)")