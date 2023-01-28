%% Settings
clearvars
close all
clc
format short
%% System Setup
s = tf('s');
fall = 'SISO'; % SISO or MIMO
switch fall
    case 'SISO'
        %% SISO
        n = 3;
        A = [3,4,-3;-6,7,1;2,1,4];
        A = [-2,1,-3;0,-3,-1;0,0,-1];
        % A = [2,-1,0;-1,-1,3;0,3,-4];
        Bsiso = [1;0;2];
        Csiso = [1,-1,2];
        Dsiso = 0;

        Gspectral = spectralDecompSSR(s,A,Bsiso,Csiso,Dsiso)
        SSR = ss(A,Bsiso,Csiso,Dsiso);
        G = minreal(tf(SSR))

    case 'MIMO'
        %% MIMO
        n = 3;
        A = [3,4,-3;-6,7,1;2,1,4];
        A = [2,1,-3;0,3,-1;0,0,-1];
        A = [2,-1,0;-1,-1,3;0,3,-4];
        Bmimo = [-1,0;4,-2;0,1];
        Cmimo = [2,0,0;1,-1,2];
        Dmimo = zeros(2,2);

        Gspectral = spectralDecompSSR(s,A,Bmimo,Cmimo,Dmimo)
        SSR = ss(A,Bmimo,Cmimo,Dmimo);
        G = minreal(tf(SSR))
end


%% Comparing responses
figure
hold on
step(Gspectral)
step(G)
grid on
