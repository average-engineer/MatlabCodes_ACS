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
n = size(A,1); % number of states
m = size(B,2); % number of inputs
l = size(C,1); % number of outputs
%% LQR Weighting Matrices
% Design 1 (Equal weightage to actuator effort and response decay)
R1 = eye(m,m); % Control Cost Matrix
Q1 = eye(n,n); % State-weighing matrix
% Design 2 (Actuator Effort is 100 times more important than response decay)
R2 = eye(m,m);
Q2 = 0.01*eye(n,n);

%% Pre-requisite for LQR
if rank(ctrb(A,B)) == n && rank(obsv(A,C)) == n
    fprintf("Pre-requisites for LQR met\n");
end
%% Solving Algebraic Riccati Equation for Optimal K
[K01,P1,clEig1] = lqr(A,B,Q1,R1)
[K02,P2,clEig2] = lqr(A,B,Q2,R2)
%% Time-Domain Response
% Initial Condition
x0 = [0;0.1;zeros(4,1)];
% Optimal Closed-Loop systems
Acl1 = A - B*K01;
Ccl1 = C - D*K01;
Acl2 = A - B*K02;
Ccl2 = C - D*K02;
sysCL1 = ss(Acl1,zeros(size(B)),Ccl1,zeros(size(D)));
sysCL2 = ss(Acl2,zeros(size(B)),Ccl2,zeros(size(D)));
[y1,t1,x1] = initial(sysCL1,x0);
[y2,t2,x2] = initial(sysCL2,x0);
for ii = 1:length(t1)
    u1(ii,:) = -K01*x1(ii,:)'; % Actuator Effort
end
for ii = 1:length(t2)
    u2(ii,:) = -K02*x2(ii,:)'; % Actuator Effort
end


figure
subplot(2,2,1)
plot(t1,y1(:,1),'DisplayName','Q = R = I','LineWidth',2)
hold on
plot(t2,y2(:,1),'DisplayName','Q = 0.01I, R = I','LineWidth',2)
legend
grid on
title("y_1(t)")

subplot(2,2,2)
plot(t1,y1(:,2),'DisplayName','Q = R = I','LineWidth',2)
hold on
plot(t2,y2(:,2),'DisplayName','Q = 0.01I, R = I','LineWidth',2)
legend
grid on
title("y_2(t)")

subplot(2,2,3)
plot(t1,u1(:,1),'DisplayName','Q = R = I','LineWidth',2)
hold on
plot(t2,u2(:,1),'DisplayName','Q = 0.01I, R = I','LineWidth',2)
legend
grid on
title("u_1(t)")

subplot(2,2,4)
plot(t1,u1(:,2),'DisplayName','Q = R = I','LineWidth',2)
hold on
plot(t2,u2(:,2),'DisplayName','Q = 0.01I, R = I','LineWidth',2)
legend
grid on
title("u_2(t)")