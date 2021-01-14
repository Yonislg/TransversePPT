%function [I_t, V_t] = current(V_0, R, L , C)
%close all
clear all
L = 38.87 * 10^-9;       %58.89
R = 24.61 / 1000;   %11.45
C = 20 * 10^-6;
V_0 = 1300;

a   = -R/(2*L);
w_0 = 1/sqrt(L*C);
dmp = a/w_0;
wd  = sqrt(1/(L*C) - R^2/(4*L^2));
    
    
I_0 = V_0/(wd*L)

T = 2*pi/wd

tmax = T/4

I_max = I_0*exp(a*tmax)

t = linspace(0,4*T, 500);

I = I_0*exp(a.*t).*sin(wd.*t);

plot(t,I)

%end