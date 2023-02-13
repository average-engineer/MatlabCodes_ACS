%% Settings
clearvars
close all
clc
format short
%% Plant Model
s = tf('s');
G = 1/(10*s + 1);

%% Paramters for Mu-Toolbox
systemnames = 'G'; % Plant Model
inputvar = '[d;r;n;u]'; % Exogeneous Inputs and Control Signals
input_to_G = '[u]'; % Controller Signal
outputvar = '[G + d + r ; r - G - d - n]'; % Exogeneous and Sensed Outputs (as a function of all inputs)
sysoutname = 'P';
P = sysic;

%% Generating a lower LFT with a controller
kp = 10;
kd = 100;
ki = 1;
K = kp + kd*s; % Random PD Controller
lowLFT = starp(evalfr(P,0),evalfr(K,0));