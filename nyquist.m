%% Settings
clearvars
close all
clc
format short
%% Plant
s = tf('s');
G = (2*s*s + 5*s + 1)/(s*s - 2*s + 3);
%% Controller
K = 1;
%% Open Loop Gain
L = G*K;
%% Nyquist Plot
figure
nyquist(L)
grid on

figure
nyquist(1 + L)
grid on