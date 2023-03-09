%% Settings
clearvars
close all
clc
format short
%% Plant
[A,B,C,D] = tf2ss(3*[-2,1],conv([5,1],[10,1]));