%% Author: Ashutosh Mukherjee
% Co-Author: Raj Khamkar
% Last Date of modification: 26.02.2023
%% Settings
clearvars
close all
clc
format short
addpath("matlabtikz")
%% System Setup
s = tf('s');
fall = 'SISO'; % SISO or MIMO
switch fall
    case 'SISO'
        %% SISO
        n = 3;
        A = [-2,1,-0.5;1,-1,0;-0.5,0,-4];
        Bsiso = [1;0;2];
        Csiso = [1,-1,2];
        Dsiso = 0;

        Gspectral = spectralDecompSSR(s,A,Bsiso,Csiso,Dsiso)
        SSR = ss(A,Bsiso,Csiso,Dsiso);
        G = minreal(tf(SSR))

        %% Comparing responses
        t = linspace(0,5,100);
        u = 2.*exp(-t);
        figure
        hold on
        x1 = lsim(Gspectral,u,t);
        x2 = lsim(G,u,t);
        plot(t,x1,'linewidth',2,'color','r','DisplayName','Spectral Decomposition')
        plot(t,x2,'--','linewidth',3,'color','k','DisplayName','MATLAB ss2tf')
        xlabel('Time (s)')
        ylabel('System Response')
        grid on
        legend
        matlab2tikz()

    case 'MIMO'
        %% MIMO
        n = 3;
        A = [-2,1,-0.5;1,-1,0;-0.5,0,-4];
        Bmimo = [-1,0;4,-2;0,1];
        Cmimo = [2,0,0;1,-1,2];
        Dmimo = zeros(2,2);

        Gspectral = spectralDecompSSR(s,A,Bmimo,Cmimo,Dmimo)
        SSR = ss(A,Bmimo,Cmimo,Dmimo);
        G = minreal(tf(SSR))

        %% Comparing responses
        t = linspace(0,5,100);
        u = 2.*exp(-t);
        uip = [u;zeros(size(t))];
        figure
        x1 = lsim(Gspectral,uip,t);
        x2 = lsim(G,uip,t);
        subplot(2,1,1)
        hold on
        plot(t,x1(:,1),'linewidth',2,'color','r','DisplayName','Spectral Decomposition y1')
        plot(t,x2(:,1),'--','linewidth',3,'color','k','DisplayName','MATLAB ss2tf y1')
        xlabel('Time (s)')
        ylabel('System Response')
        grid on
        legend
        subplot(2,1,2)
        hold on
        plot(t,x1(:,2),'linewidth',2,'color','r','DisplayName','Spectral Decomposition y2')
        plot(t,x2(:,2),'--','linewidth',3,'color','k','DisplayName','MATLAB ss2tf y2')
        xlabel('Time (s)')
        ylabel('System Response')
        grid on
        legend
        matlab2tikz()
end



