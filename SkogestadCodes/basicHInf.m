clearvars
close all
clc
format short
% Uses the Mu-toolbox
G=nd2sys(1,conv([10 1],conv([0.05 1],[0.05 1])),200); % Plant is G.
M=1.5; wb=10; A=1.e-4; Wp = nd2sys([1/M wb], [1 wb*A]); Wu = 1; % Weights.
systemnames = 'G Wp Wu';
inputvar = '[ r(1); u(1)]';
outputvar = '[Wp; Wu; r-G]';
input_to_G = '[u]';
input_to_Wp = '[r-G]';
input_to_Wu = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
%
% Find H-infinity optimal controller:
%
nmeas=1; nu=1; gmn=0.5; gmx=20; tol=0.001;
[khinf,ghinf,gopt] = hinfsyn(P,nmeas,nu,gmn,gmx,tol);