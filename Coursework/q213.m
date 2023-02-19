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
part = 'c';
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

        %% Upper Singular Value Frequency Response
        for kk = 1:length(freq)
            % Additive Uncertainty Singular Value
            [~,eeA,~] = svd(evalfr(MA,freq(kk)*1i));
            MAFR(kk) = eeA(1,1);
            % Multiplicative Uncertainty Singular Value
            [~,eeM,~] = svd(evalfr(MM,freq(kk)*1i));
            MMFR(kk) = eeM(1,1);
        end

        % Computing the frequency at which Peak is reached
        [HInfNormA,peakFreqA] = hinfnorm(MA);
        [HInfNormM,peakFreqM] = hinfnorm(MM);

        figure
        semilogx(freq,MAFR,'-*','DisplayName','Additive Uncertainty')
        hold on
        semilogx(freq,MMFR,'-*','DisplayName','Output Multiplicative Uncertainty')
        grid on
        xlabel('Frequency (rad/s)')
        title('Singular Value Plot')
        legend 
        display(append("H-Inf Norm of MM is ",num2str(hinfnorm(MM))));
        display(append("H-Inf Norm of MA is ",num2str(hinfnorm(MA))));

    case 'c'
        MA = minreal(K*inv(eye(size(G)) - G*K));
        MM = minreal(G*K*inv(eye(size(G)) - G*K));
        gamma = 0.06; % Worst Case Perturbation
        critFreq = 2.83; % Critical Frequency (rad/s) at which worst case perturbation occurs
        %% Lower Singular Value of Perturbation Matrix at critical frequency
        for kk = 1:length(freq)
            [~,eeA,~] = svd(evalfr(MA,freq(kk)*1i));
            lowMA(kk) = eeA(end,end);
            [~,eeM,~] = svd(evalfr(MM,freq(kk)*1i));
            lowMM(kk) = eeM(end,end);
        end
        %% Perturbation Matrix

        

       
        
end


