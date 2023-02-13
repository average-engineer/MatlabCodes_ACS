%% Settings
clearvars
close all
clc
format short
s = tf('s');
%% Plant
G = [7,8;6,7]*[1/(s+1),0;0,2/(s+2)]*inv([7,8;6,7]);
%% Controller
K = -eye(size(G));
LI = K*G;
L = G*K;
olTF = eye(size(G)) - L;
detTF = (olTF(1,1)*olTF(2,2) - olTF(1,2)*olTF(2,1));
%% Closed Loop TFs
Tf1 = minreal(inv(eye(size(G)) - (K*G)));
Tf2 = minreal(K*inv(eye(size(G)) - (G*K)));
Tf3 = minreal(G*Tf1);
Tf4 = minreal(inv(eye(size(G)) - (G*K)));
%% Question Part
part = 'a';
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
        freq = logspace(-3,2,100);
        %% Lower LFT (Manual)
        P11 = [zeros(size(G)),eye(size(G)),zeros(size(G))];
        P12 = eye(size(G));
        P21 = [eye(size(G)),G,eye(size(G))];
        P22 = G;
        P = [P11,P12;P21,P22];
        N = lft(P,K);
        M = N(1:2,1:2);
        for kk = 1:length(freq)
            Mf = evalfr(M,freq(kk)*1i);
            [~,eeM,~] = svd(Mf);
            MFR(kk) = eeM(1,1);
        end
        
        %% Perturbatons
        del = linspace(-0.1,2,50);
        % Additive Uncertainty
        countA = 1;
        for jj = 1:length(del)
            % delM = del(jj)*eye(size(G));
            delM = [del(jj),-del(jj);-del(jj),del(jj)];
            LpA(:,:,countA) = L + delM;
            LpM(:,:,countA) = L*(eye(size(G)) + delM);
            for kk = 1:length(freq)
                LpfA = evalfr(LpA(:,:,countA),freq(kk)*1i);
                [~,eeLpA,~] = svd(LpfA);
                LpAFR(kk,countA) = eeLpA(1,1);
            end
            if max(LpAFR(:,countA)) >= 1
                % countA = countA-1;
                break;
            end
            countA = countA + 1;

        end

        % Multiplicative Output Uncertainty
        countM = 1;
        for jj = 1:length(del)
            % delM = del(jj)*eye(size(G));
            delM = [del(jj),-del(jj);-del(jj),del(jj)];
            LpA(:,:,countM) = L + delM;
            LpM(:,:,countM) = L*(eye(size(G)) + delM);
            for kk = 1:length(freq)
                LpfM = evalfr(LpM(:,:,countM),freq(kk)*1i);
                [~,eeLpM,~] = svd(LpfM);
                LpMFR(kk,countM) = eeLpM(1,1);
            end
            if max(LpMFR(:,countM)) >= 1
                % countM = countM-1;
                break;
            end
            countM = countM + 1;

        end
        
        
        %% Frequency Response 
        for kk = 1:length(freq)
            Lf = evalfr(L,freq(kk)*1i);
            [~,eeL,~] = svd(Lf);
            LFR(kk) = eeL(1,1);
        end

        figure
        semilogx(freq,LFR,'--','color','k','linewidth',2,'DisplayName','Open-Loop Singular Value')
        hold on
        for jj = 1:countA
            semilogx(freq,LpAFR(:,jj),'DisplayName',append('Perturbation = ',num2str(del(jj))))
        end
        grid on
        title('Additive Uncertainty')
        legend

        figure
        semilogx(freq,LFR,'--','color','k','linewidth',2,'DisplayName','Open-Loop Singular Value')
        hold on
        for jj = 1:countM
            semilogx(freq,LpMFR(:,jj),'DisplayName',append('Perturbation = ',num2str(del(jj))))
        end
        grid on
        title('Multiplicative Output Uncertainty')
        legend

         
        
        
end


