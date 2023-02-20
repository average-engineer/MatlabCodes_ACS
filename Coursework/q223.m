%% Settings
clearvars
close all
clc
format short
%% Verification of Bezout Identity
N = tf(zpk(1,[-1 -1 -1],1));
M = tf(zpk([-2 2],[-1 -1],1));
U = tf([61 121],[3 3]);
V = tf([3 12 -31],[3 6 3]);
bzIden = U*N + V*M
%% Youla Stabilizing Controller
Q = tf(1,[1 3]);
K = (U + Q*M)/(V - Q*N)