%% Settings
clearvars
close all
clc
format short

%% LTI System (TISO System)
A = [8,3;-1,4];
B = [1,0;0,2];
C = [1,2];
D = 0;

sys = ss(A,B,C,D);

%% Conanical Forms
% Expected A in:
% Jordon Form: [5,0;0,7];
% Controller Companion: [0,1;-35,12]
% Observer Companion: [0,-35;1,12]

jordon = canon(sys,'model');

obs = canon(sys,'companion');