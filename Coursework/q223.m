%% Settings
clearvars
close all
clc
format short
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
%% Step Response
figure
subplot(2,2,1)
step(Tf1)
title('Tf1')
grid on
subplot(2,2,2)
step(Tf2)
title('Tf2')
grid on
subplot(2,2,3)
step(Tf3)
title('Tf3')
grid on
subplot(2,2,4)
step(Tf4)
title('Tf4')
grid on