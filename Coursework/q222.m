%% Settings
close all
clearvars
clc
format short
s = tf('s');
%% State Space Realisation (Minimal)
A = [2,0;0,-2];
B = eye(2,2);
C = [1,-1];
D = [1,1];
SYS = ss(A,B,C,D);
%% Polynomial Matrix
P = [s*eye(size(A)) - A,-B;C,D];
%% Zeros
z = tzero(SYS);