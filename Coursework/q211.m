%% Settings
clearvars
close all
clc
format short
% s = tf('s');
s = 5; % Value of s at which validation is done
omega = linspace(0,100);
% s = 1i*omega;
%% Inverted Pendulum on Cart
n = 4; % Number os states
% States
% x1: Cart Horizontal Position
% x2: Pendulum Angle with the vertical
% x3: Cart Horizontal Speed
% x4: Pendulum angular veloctiy
mCart = 5; % Mass of Cart
mPend = 1; % Mass of Pendulum
L = 1; % Length of Pendulum
g = 9.81; % Gravity
% Governing Linearised Equations of Motion
M = [mCart + mPend,mPend*L/2;mPend*L/2,mPend*L*L/3]; % Mass Matrix
P = zeros(size(M)); % Damping Matrix
Q = [0,0;0,-mPend*g*L/2]; % Stiffness Matrix
A = [zeros(n/2,n/2),eye(n/2,n/2);-M\Q,-M\P];
%% SISO Realisation
% % System Input -> Motor actuation on wheels of cart
% % System Output -> Angle made by the pendulum along the vertical
% F = 10; % Wheel Motor Actuation
% h = [F;0]; % Excitation Vector
% B = [zeros(n/2,1);M\h];
% Bsiso = B/F;
% Csiso = [0,1,0,0];
% Dsiso = 0;
% 
% Gspectral = spectralDecompSSR(s,A,Bsiso,Csiso,Dsiso);
% 
% % Transfer Function Obtained from taking laplace transform of the EOMs
% % Proof in report
% % GLaplace = 6/(6*(mCart + mPend)*g - L*s*s*(4*mCart + mPend));
% 
% G = Csiso*((s*eye(n,n) - A)\Bsiso) + Dsiso;

%% Random Plant (SISO)
n = 3;
A = [3,4,-3;-6,7,1;2,1,4];
Bsiso = [1;0;2];
Csiso = [1,-1,2];
Dsiso = 0;

% for k = 1:length(omega)
%     Gspectral(k) = spectralDecompSSR(s(k),A,Bsiso,Csiso,Dsiso);
%     G(k) = Csiso*((s(k)*eye(n,n) - A)\Bsiso) + Dsiso;
% end

% figure
% hold on
% plot(omega,abs(Gspectral),'linewidth',2)
% plot(omega,abs(G),'--','linewidth',2,'color','k')
% grid on

%% Random Plant (MIMO)
n = 3;
A = [3,4,-3;-6,7,1;2,1,4];
Bmimo = [-1,0;4,-2;0,1];
Cmimo = [2,0,0;1,-1,2];
Dmimo = zeros(2,2);

Gspectral = spectralDecompSSR(s,A,Bmimo,Cmimo,Dmimo);
G = Cmimo*((s*eye(n,n) - A)\Bmimo) + Dmimo;