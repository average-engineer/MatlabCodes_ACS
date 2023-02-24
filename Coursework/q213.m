%% Settings
clearvars
close all
clc
format short
s = tf('s');
freq = logspace(-3,3,500);
%% Plant
G = [7,8;6,7]*[1/(s+1),0;0,2/(s+2)]*inv([7,8;6,7]);
%% Controller
K = -eye(size(G));
LI = K*G;
L = G*K;
%% Closed Loop TFs
Tf1 = minreal(inv(eye(size(G)) - (K*G)));
Tf2 = minreal(K*inv(eye(size(G)) - (G*K)));
Tf3 = minreal(G*Tf1);
Tf4 = minreal(inv(eye(size(G)) - (G*K)));
%% Question Part
part = 'b'; % a or b
switch part
    case 'a'
        %% Pole Zero Map
        figure
        subplot(2,2,1)
        pzmap(Tf1)
        grid on
        subplot(2,2,2)
        pzmap(Tf2)
        grid on
        subplot(2,2,3)
        pzmap(Tf3)
        grid on
        subplot(2,2,4)
        pzmap(Tf4)
        grid on
        poles = pole(Tf1)
        %% Arbritary Input
        t = linspace(0,10,100);
        u1 = 2.*exp(-t);
        u2 = zeros(size(t));
        u = [u1;u2];
        figure
        subplot(2,2,1)
        y1 = lsimplot(Tf1,u,t);
        grid on
        title('Tf1')
        subplot(2,2,2)
        y2 = lsimplot(Tf2,u,t);
        grid on
        title('Tf2')
        subplot(2,2,3)
        y3 = lsimplot(Tf3,u,t);
        grid on
        title('Tf3')
        subplot(2,2,4)
        y4 = lsimplot(Tf4,u,t);
        grid on
        title('Tf4')

        figure
        lsimplot(G,u,t)
        grid on
        title('Plant')
        
        freq = logspace(-3,2,100);
        for kk = 1:length(freq)
            Lf = evalfr(L,freq(kk)*1i);
            [~,eeL,~] = svd(Lf);
            LFR(kk) = eeL(1,1);
        end
        figure
        semilogx(freq,LFR,'--','color','k','linewidth',2,'DisplayName','Open-Loop Singular Value')
        grid on
        legend

    case 'b'
        
        %% Additive Dynamic Uncertainty
        % M Matrix
        MA = minreal(K*inv(eye(size(G)) - G*K));

        %% Multiplicative Output Dynamic Uncertainty
        % M Matrix
        MM = minreal(G*K*inv(eye(size(G)) - G*K));

        %% Combined Dynamic Uncertainty
        % Generalized Plant assuming combined additive and op
        % multiplicative uncertainty
        P = [zeros(2,2),zeros(2,2),eye(2,2),zeros(2,2),eye(2,2);
            eye(2,2),zeros(2,2),G,zeros(2,2),G;
            eye(2,2),eye(2,2),G,zeros(2,2),G;
            eye(2,2),eye(2,2),G,eye(2,2),G];
        N = lft(P,K);
        M = N(1:4,1:4);

        %% Upper Singular Value Frequency Response
        for kk = 1:length(freq)
            [~,ee,~] = svd(evalfr(M,freq(kk)*1i));
            MFR(kk) = ee(1,1);
        end

        % Computing the frequency at which Peak is reached
        [HInfNorm,peakFreq] = hinfnorm(M);

        figure
        semilogx(freq,MFR,'-*','DisplayName','Dynamic Uncertainty M Upper Singular Value')
        hold on
        grid on
        xlabel('Frequency (rad/s)')
        title('Singular Value Plot')
        legend 
        fprintf(append("H-Inf Norm of M is ",num2str(HInfNorm)," and occurs at ",num2str(peakFreq), " rad/s\n"));

        gamma = round(HInfNorm);
        %% Perturbation Matrix
        delA = gamma*tf([1 -2],[1 2]);
        delM = 0.01*tf([1 -2],[1 2]);
        pertM = [delA*eye(2,2),zeros(2,2);zeros(2,2),delM*eye(2,2)];
        for kk = 1:length(freq)
            [~,eeP,~] = svd(evalfr(pertM,freq(kk)*1i));
            PertFR(kk) = eeP(1,1);
        end

        figure
        semilogx(freq,PertFR,'-*','DisplayName','Perturbation Matrix Upper Singular Value')
        hold on
        grid on
        xlabel('Frequency (rad/s)')
        title('Singular Value Plot')
        legend 


        %% Closing the upper loop with perturbation matrix
        F = lft(pertM,N);
        % Step Response 
        figure
        step(F)
        grid on
        title('Interconnection of closed-loop and perturbation')      
        
end


