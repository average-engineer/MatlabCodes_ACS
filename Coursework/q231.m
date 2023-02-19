%% Settings
clearvars
close all
clc
format short
s = tf('s');
freq = logspace(-3,3,100);
fall = 1;
%% Verifying the model for parametric uncertainties
mnom = 10;
mmax = 18;
mmin = 5;
bm = 0.385*mnom;
dm = -0.23;
delm = linspace(-1,1,21);
for ii = 1:length(delm)
    m(ii) = (mnom + bm*delm(ii))/(1 + dm*delm(ii));
end

knom = 1;
kmax = 1.8;
kmin = 0.8;
bk = -0.28*knom;
dk = -0.6;
delk = linspace(-1,1,21);
for ii = 1:length(delk)
    k(ii) = (knom + bk*delk(ii))/(1 + dk*delk(ii));
end

hnom = 2.5;
bh = -0.28*hnom;
dh = -0.6;
cnom = 10;
cmax = 14;
cmin = 6;
rc = (cmax-cmin)/(cmax + cmin);

%% Defining Parametric Uncertainties
mp = ureal('m',mnom,'Percentage',[-50,80]);
cp = ureal('c',cnom,'Percentage',[-40,40]);
kp = ureal('k',knom,'Percentage',[-20,80]);
hp = ureal('h',hnom,'Percentage',[-20,80]);
%% Weights
W1 = tf([2 0],[1 0.1127]); % Output Weight Weight
W2 = tf([0.7 1],[0.07 1]); % Actuator Effort Weight
W2 = 0.3*W2;
%% PI Controller Design for Nominal System
K0 = 10*tf([10 1],[10 0]); % Lag Compensator
K = W2*K0; % W2 acting as a Lead Compensator
Lnom = minreal(W1*hnom*K*tf(1,[mnom cnom knom]));
switch fall
    case 1

        %% Uncertain Plant Transfer Function
        % Input: F
        % Output: y
        uPlant = tf(hp,[mp cp kp]);
        uPlant = W1*uPlant;
        %% Nominal Stability
        clTF1 = W1/(1 + Lnom); % closed loop TF from disturbance to output
        for kk = 1:length(freq)
            LnomFR(kk) = abs(evalfr(Lnom,freq(kk)*1i));
        end

        %         figure
        %         semilogx(freq,LnomFR,'--','color','k','DisplayName','Lead-Lag Compensated Loop Gain')
        %         xlabel('Frequency (rad/s)')
        %         grid on
        %         legend

        figure
        step(clTF1)
        grid on

        %% Manually Obtained Generalized Plant
        P = [-bm/mnom,-1,-(bk - knom*dk),0,-knom,-cnom,0,1;
            zeros(1,5),rc*cnom,0,0;
            0,0,-dk,0,1,zeros(1,3);
            zeros(1,3),-dh,1,zeros(1,3);
            zeros(1,5),1,0,0;
            (dm/mnom) - (bm/(mnom^2)),-1/mnom,-(bk-knom*dk)/mnom,0,-knom/mnom,-cnom/mnom,0,1/mnom;
            zeros(1,3),W1*(bh-hnom*dh),W1*hnom,0,W1,0;
            zeros(1,3),-W1*(bh-hnom*dh),-W1*hnom,0,-W1,0];
        %% Closing the Lower Controller Loop
        P11 = P(1:7,1:7);
        P12 = P(1:7,8);
        P21 = P(8,1:7);
        P22 = P(8,8);
        % N = P11 + P12*K*inv(eye(size(P22*K)) - P22*K)*P21;
        N = lft(P,K);
        M = minreal(N(1:4,1:4));

        %% Using Robust Control Toolbox
        % Uncertain Closed Loop
        Lp = W1*hp*K*tf(1,[mp cp kp]);
        clTFp = uss(W1/(1 + Lp));
        % Uncertain State Space (Open Loop)
        Aol = [0,1;-kp/mp,-cp/mp];
        Bol = [0;1/mp];
        Col = [hp,0];
        Dol = 0;
        % Uncertain State Space (Closed Loop)
        Acl = Aol - Bol*K*W1*[hp,0];
        Bcl = -Bol*K*W1;
        Ccl = [W1*hp,0];
        Dcl = W1;
        CLsys = ss(Acl,Bcl,Ccl,Dcl);
        [N1,delta] = lftdata(CLsys);
        N1 = minreal(tf(N1));
        M1 = N1(1:size(delta,1),1:size(delta,1));

        %% Robust Stability
        % Perturbation Block Structure
        % All 4 Perturbation Blocks: Real single perturbation
        BlockStructure = [-1,0;-1,0;-1,0;-1,0];
        for kk = 1:length(freq)
            % Using upper singular value
            [~,ee,~] = svd(evalfr(M,freq(kk)*1i));
            [~,ee1,~] = svd(evalfr(M1,freq(kk)*1i));
            MSVFR(kk) = ee(1,1);
            M1SVFR(kk) = ee1(1,1);
            % Using Structure Singular Value
            bounds = mussv(evalfr(M,freq(kk)*1i),BlockStructure);
            MSSVFR(kk) = max(bounds); % Using the upper bound (more reliable due to convex optimization)
            bounds1 = mussv(evalfr(M1,freq(kk)*1i),BlockStructure);
            M1SSVFR(kk) = max(bounds1); % Using the upper bound (more reliable due to convex optimization)
        end

        figure
        semilogx(freq,MSVFR,'--','color','k','linewidth',2,'DisplayName','Upper Singular Value (M)')
        hold on
        semilogx(freq,M1SVFR,'color','r','DisplayName','Upper Singular Value (M1)')
        xlabel('Frequency (rad/s)')
        grid on
        legend
        figure
        semilogx(freq,MSSVFR,'--','color','k','linewidth',2,'DisplayName','Structured Singular Value (M)')
        hold on
        semilogx(freq,M1SSVFR,'color','r','DisplayName','Structured Singular Value (M1)')
        xlabel('Frequency (rad/s)')
        grid on
        legend

        %% Nominal Performance
        N22 = minreal(N(end-2:end,end-2:end));
        for kk = 1:length(freq)
            [~,ee,~] = svd(evalfr(N22,freq(kk)*1i));
            N22FR(kk) = ee(1,1);
        end
        figure
        semilogx(freq,N22FR,'--','color','k','DisplayName','Upper Singular Value (N22)')
        xlabel('Frequency (rad/s)')
        grid on
        legend

    case 2
        %% Including Uncertain Time Delay (Dynamic Uncertainty)
        taunom = 0.1;
        taumax = 1.3*taunom;
        % Uncertainyty weight
        wTau = tf([2*taumax 0],[taumax 2]);
        % Complex Perturbation for time delay
        deltau = ultidyn('delT',[1,1]);
        %% Uncertain Plant Transfer Function
        % Input: F
        % Output: y
        uPlant = W2*tf(hp,[mp cp kp]);
        uPlant = uPlant*(1 + wTau*deltau);
        %% Manually Obtained Generalized Plant
        P = [-bm/mnom,-1,-(bk - knom*dk),0,1,-knom,-cnom,0,W2;zeros(1,6),rc*cnom,0,0;...
            0,0,-dk,0,0,1,zeros(1,3);zeros(1,3),-dh,0,1,zeros(1,3);zeros(1,8),wTau*W2;zeros(1,6),1,0,0;...
            (dm/mnom) - (bm/(mnom^2)),-1/mnom,-(bk-knom*dk)/mnom,0,1/mnom,-knom/mnom,-cnom/mnom,0,W2/mnom;...
            zeros(1,3),W1*(bh-hnom*dh),0,W1*hnom,0,W1,0;zeros(1,3),-W1*(bh-hnom*dh),0,-W1*hnom,0,-W1,0];
        %% Closing the Lower Controller Loop
        N = lft(P,K);
        P11 = P(1:8,1:8);
        P12 = P(1:8,9);
        P21 = P(9,1:8);
        P22 = P(9,9);
        M = minreal(N(1:5,1:5));

        %% Robust Stability
        % Perturbation Block Structure
        % First 4 Perturbation Blocks: Real scalar perturbation
        % Last Perturbation Block: Complex scalar perturbation
        BlockStructure = [-1,0;-1,0;-1,0;-1,0;1,0];
        for kk = 1:length(freq)
            % Using upper singular value
            [~,ee,~] = svd(evalfr(M,freq(kk)*1i));
            MSVFR(kk) = ee(1,1);
            % Using Structure Singular Value
            bounds = mussv(evalfr(M,freq(kk)*1i),BlockStructure);
            MSSVFR(kk) = max(bounds); % Using the upper bound (more reliable due to convex optimization)
        end

        figure
        semilogx(freq,MSVFR,'--','color','k','DisplayName','Upper Singular Value (M)')
        xlabel('Frequency (rad/s)')
        grid on
        legend
        figure
        semilogx(freq,MSSVFR,'--','color','r','DisplayName','Structured Singular Value (M)')
        xlabel('Frequency (rad/s)')
        grid on
        legend

        %% Using Robust Control Toolbox
        % Uncertain Closed Loop
        Lp = W1*hp*W2*K*tf(1,[mp cp kp])*(1 + wTau*deltau);
        clTFp = uss(W1/(1 + Lp));
        opts = robOptions('VaryFrequency','on');
        [margins,wc,info] = robstab(clTFp,{freq(1),freq(end)},opts);
        figure
        semilogx(info.Frequency,info.Bounds)
        title('Stability Margin vs. Frequency')
        ylabel('Margin')
        xlabel('Frequency')
        legend('Lower bound','Upper bound')

        %% Nominal Performance
        N22 = minreal(N(end-2:end,end-2:end));
        for kk = 1:length(freq)
            [~,ee,~] = svd(evalfr(N22,freq(kk)*1i));
            N22FR(kk) = ee(1,1);
        end
        figure
        semilogx(freq,N22FR,'--','color','k','DisplayName','Upper Singular Value (N22)')
        xlabel('Frequency (rad/s)')
        grid on
        legend



end




