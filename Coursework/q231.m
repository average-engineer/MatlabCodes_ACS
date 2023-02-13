%% Settings
clearvars
close all
clc
format short
s = tf('s');
%% Parametric Uncertainities
mnom = 10;
mmax = 18;
mmin = 5;
rm = (mmax-mmin)/(mmax+mmin);
delm = linspace(-1,1,21);
for ii = 1:length(delm)
    m(ii) = mnom*(1 + rm*delm(ii));
end

knom = 1;
kmax = 1.8;
kmin = 0.8;
rk = (kmax-kmin)/(kmax+kmin);
delk = linspace(-1,1,21);
for ii = 1:length(delk)
    k(ii) = knom/(1 - rk*delk(ii));
end
