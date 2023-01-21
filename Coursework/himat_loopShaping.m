% Matlab Himat loop-shaping procedure

% Taken from Matlab demos
% Author: B. Misgeld
% Date:   2012-05-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ag =[
-2.2567e-02  -3.6617e+01  -1.8897e+01  -3.2090e+01   3.2509e+00  -7.6257e-01;
9.2572e-05  -1.8997e+00   9.8312e-01  -7.2562e-04  -1.7080e-01  -4.9652e-03;
1.2338e-02   1.1720e+01  -2.6316e+00   8.7582e-04  -3.1604e+01   2.2396e+01;
0            0   1.0000e+00            0            0            0;
0            0            0            0  -3.0000e+01            0;
0            0            0            0            0  -3.0000e+01];
bg = [0     0;
      0     0;
      0     0;
      0     0;
     30     0;
      0    30];
cg = [0     1     0     0     0     0;
      0     0     0     1     0     0];
dg = [0     0;
      0     0];
G=ss(ag,bg,cg,dg);

%% Enter the desired loop shape Gd
s=zpk('s'); % Laplace variable s
Gd=8/s; % desired loop shape
%% Compute the optimal loop shaping controller K
[K,CL,GAM]=loopsyn(G,Gd);
%% Compute the loop L, sensitivity S and 
%% complementary sensitivity T:
L=G*K;
I=eye(size(L));
S=feedback(I,L); % S=inv(I+L);
T=I-S;
%% Plot the results:
% step response plots
step(T);title('\alpha and \theta command step responses');
% frequency response plots
figure;
sigma(I+L,'--',T,':',L,'r--',Gd,'k-.',Gd/GAM,'k:',...
 Gd*GAM,'k:',{.1,100});grid
legend('1/\sigma(S) performance',...
 '\sigma(T) robustness',...
 '\sigma(L) loopshape',...
 '\sigma(Gd) desired loop',...
 '\sigma(Gd) \pm GAM, dB');