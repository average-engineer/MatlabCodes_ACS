%% Settings
clearvars
close all
clc
format short
s = tf('s');
tol = 0.1; % tolerance for difference between Relative error and weight
%% Nominal Plant
Gnom = 1/(s+1);
%% Multiplicative Uncertainty Weight
wM = (0.125*s + 0.25)/((0.125/4)*s + 1);

%% Uncertain Plant Families
fall = 'f';
switch fall
    case 'a'
        parameter = linspace(0.1,0.5,20);
        Plant = Gnom*exp(-parameter*s);

    case 'b'
        parmeter = linspace(0.1,0.5,20);
        for kk = 1:length(parameter)
            Plant(kk) = Gnom/(parameter(kk)*s + 1);
        end

    case 'c'
        parameter = linspace(0.5,1.5,20);
        for kk = 1:length(parameter)
            Plant(kk) = 1/(s+parameter(kk));
        end

    case 'd'
        parameter = linspace(0.5,1.5,20);
        for kk = 1:length(parameter)
            Plant(kk) = 1/(parameter(kk)*s+1);
        end
    case 'f'
        parameter = linspace(1,20,20);
        for kk = 1:length(parameter)
            Plant(kk) = Gnom*((1/(0.01*s+1))^parameter(kk));
        end
end
relError = (Plant - Gnom)/Gnom;

%% Comparing Frequency Repsonse
freq = logspace(-2,2,100);

for kk = 1:length(freq)
    wMFR(kk) = abs(evalfr(wM,freq(kk)*1i)); % Frequency response of the uncertainty weight
    relErrFR(kk,:) = abs(evalfr(relError,freq(kk)*1i)); % Frequency response of the Uncertain Plant
end

figure
semilogx(freq,20*log10(wMFR),'--','linewidth',2,'color','k','DisplayName','Uncertainty Weight')
hold on
for ii = 1:size(Plant,2)
    semilogx(freq,20*log10(relErrFR(:,ii)),'linewidth',1,'DisplayName',append('Relative Uncertainty for parameter value = ',num2str(parameter(ii))))
end

grid on
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
legend
