%% Settings
clearvars
close all
clc
%% Plant
s = tf('s');
G = (1/(75*s + 1))*[87.8,-86.4;108.2,-109.6]; % Distillation Process

%% Controller
design = 'b';
switch design
    case 'a'
        c1 = 0.005;
        c2 = 0.005;
    case 'b'
        c1 = 0.005;
        c2 = 0.05;
    case 'c'
        c1 = 0.0036;
        c2 = 0.504;
end
Ks = [c1*((75*s + 1)/s),0;0,c2*((75*s + 1)/s)];

[U,~,V] = svd(evalfr(G,0)); % SVD of plant at steady state (s = 0)
W1 = V;
W2 = U.';
K = W1*Ks*W2; % SVD Controller

%% Closed Loop System
L = G*K; % Loop Gain
T = (eye(size(L)) + L)\L; % Complimentary Sensitivity TF
T = minreal(T);
S = inv(eye(size(L)) + L); % Sensitivity TF
S = minreal(S);

%% Closed Loop Response
t = 0:0.1:60; % Time vector
t = t';
r1 = 0.2*exp(-0.2*t); % r1 = 1/(5*s + 1), r2 = 0
r2 = zeros(size(t));
r = [r1,r2]; % Reference input

figure
lsimplot(T,r,t)
title(append('Design ',design))
grid on

%% Modelling input channel uncertainities
e1 = 0.2;
e2 = -0.2;

%% Perturbed Plant and Closed Loop System
Gper = G*[1+e1,0;0,1+e2];
Lper = minreal(Gper*K);
Tper = (eye(size(Lper)) + Lper)\Lper;
Tper = minreal(Tper);

%% Perturbed System Response
figure
lsimplot(Tper,r,t)
title(append('Design ',design))
grid on

%% Modelling output channel uncertainities
ee1 = e2;
ee2 = e1;
GperIO = [1+ee1,0;0,1+ee2]*G*[1+e1,0;0,1+e2];
LperIO = minreal(GperIO*K);
TperIO = (eye(size(LperIO)) + LperIO)\LperIO;
TperIO = minreal(TperIO);

%% Perturbed System Response
figure
lsimplot(TperIO,r,t)
title(append('Design ',design))
grid on



