%% Settings
clearvars
close all
clc
format short

%% Uncertain Parameter Ranges
% Controller Gain
kmin = 2;
kmax = 3;
% 
% Time Delay
thetaMin = 2;
thetaMax = 3;
% 
% Time Constant
tauMin = 2;
tauMax = 3;
% 
%% Plant Family
s = tf('s');

generation = 'loop'; % usample or loop

switch generation
    case 'loop'
        %% Using Loops
        k = [kmin:0.1:kmax];
        theta = [thetaMin:0.1:thetaMax];
        tau = [tauMin:0.1:tauMax];
        count = 1;
        for ii = 1:length(k)
            for jj = 1:length(theta)
                for kk = 1:length(tau)
                    Plants(count) = (k(ii)/(tau(kk)*s + 1))*exp(-theta(jj)*s);
                    count = count + 1;
                end
            end
        end
        point = 0.2*1i; % selected point
        res = evalfr(Plants,point);

    case 'usample'
        % Defining uncertain parameters
        k = ureal('k',(kmax + kmin)/2,'Range',[kmin,kmax]);
        theta = ureal('theta',(thetaMin + thetaMax)/2,'Range',[thetaMin,thetaMax]);
        tau = ureal('tau',(tauMin + tauMax)/2,'Range',[tauMin,tauMax]);

        % Implementing a 2nd Pade Approximation for time delay
        timeDelay = ((theta*theta/12)*s*s - (theta/2)*s + 1)/((theta*theta/12)*s*s + (theta/2)*s + 1);


        Gp = (k/(tau*s + 1))*timeDelay;

        Plants = usample(Gp,1000);

        point = 0.2*1i; % selected point
        res = evalfr(Plants,point);
        res = reshape(res,size(res,3),1);
end



%% Plotting the actual uncertainty region
figure
plot(res,'*')
grid on
xlabel('Re')
ylabel('Im')
title(append('Uncertainty Region at ',num2str(imag(point)),' rad/s'))

%% Nominal Plant
knom = (kmax + kmin)/2;
taunom = (tauMax + tauMin)/2;
thetanom = (thetaMax + thetaMin)/2;
Gnom = (knom/(taunom*s + 1))*exp(-thetanom*s);
%% Maximum encompassing radius
% Frequency Vector
freq = logspace(-1,1,100);
for kk = 1:length(freq)
    relError(:,kk) = abs((evalfr(Plants,freq(kk)*1i) - evalfr(Gnom,freq(kk)*1i))/evalfr(Gnom,freq(kk)*1i));
    radius(kk) = max(relError(:,kk));
end

%% Frequency Response
figure
hold on
semilogx(freq,20*log10(radius),'color','b','linewidth',4)
for ii = 1:size(relError,1)
    semilogx(freq,20*log10(relError(ii,:)),'--','color','k')
end
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Frequency Response')
grid on
legend('Relative Error Radius')
axis([freq(1),freq(end),-60,30])


%% Determining complex weight
% First Order Weight
T = 4;
weight1 = (T*s + 0.2)/((T/2.5)*s + 1);

% Second Order Weight
weight2 = weight1*((s*s + 1.6*s + 1)/(s*s + 1.4*s + 1));

% Weight obtained from FRD magnitude fitting of radius frequency response data
radiusFRD = frd(radius,freq);
weight3 = fitmagfrd(radiusFRD,2); % 2nd order weight

for kk = 1:length(freq)
    wM1(kk) = evalfr(weight1,freq(kk));
    wM2(kk) = evalfr(weight2,freq(kk));
    wM3(kk) = evalfr(weight3,freq(kk));
end

%% Comaprison of weights
figure
hold on
semilogx(freq,20*log10(radius),'--','color','b','linewidth',4)
semilogx(freq,20*log10(wM1),'color','r','linewidth',4)
semilogx(freq,20*log10(wM2),'color','g','linewidth',4)
semilogx(freq,20*log10(wM3),'color','m','linewidth',4)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Frequency Response')
grid on
legend('Relative Error Radius','First Order Weight','Second Order Weight','FRD Magntiude Fitted Weight')
axis([0,10,-60,30])

%% Disc Uncertainty visualisation
% Frequency : 0.2 rad/s
figure
hold on
% Discs due to various multiplicative weights
disc1 = nsidedpoly(1000,'Center',[real(evalfr(Gnom,point)),imag(evalfr(Gnom,point))],'Radius',abs(evalfr(Gnom*weight1,point)));
disc2 = nsidedpoly(1000,'Center',[real(evalfr(Gnom,point)),imag(evalfr(Gnom,point))],'Radius',abs(evalfr(Gnom*weight2,point)));
disc3 = nsidedpoly(1000,'Center',[real(evalfr(Gnom,point)),imag(evalfr(Gnom,point))],'Radius',abs(evalfr(Gnom*weight3,point)));
plot(disc1,'EdgeColor','r','FaceColor','r')
plot(disc2,'EdgeColor','g','FaceColor','g')
plot(disc3,'EdgeColor','m','FaceColor','m')
plot(evalfr(Plants,point),'*') % Actual Uncertainty Region
grid on
xlabel('Re')
ylabel('Im')
title(append('Uncertainty Region at ',num2str(imag(point)),' rad/s'))
legend('Disc Uncertainty (Weight 1)','Disc Uncertainty (Weight 2)','Disc Uncertainty (Weight 3)','Actual Uncertainty')

