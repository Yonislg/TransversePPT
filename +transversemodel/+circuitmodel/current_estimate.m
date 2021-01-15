%function [I_t, V_t] = current(V_0, R, L , C)
%close all
clear all
L = 38.87 ;       %58.89 * 10^-9; %
L = L* 10^-9;
R = 24.61 ;   % 11.45 /1000; %
R = R / 1000;
C = 80; % 20 ;
C = C * 10^-6;
V_0 = 1300;% sqrt(2*8/C) %1300;

a   = -R/(2*L);
w_0 = 1/sqrt(L*C);
dmp = a/w_0;
wd  = sqrt(1/(L*C) - R^2/(4*L^2));
    
    
I_0 = V_0/(wd*L)

T = 2*pi/wd

tmax = T/4

I_max = I_0*exp(a*tmax)

t = linspace(0,35*10^-6, 500);

I = I_0*exp(a.*t).*sin(wd.*t);
figure(1)
plot(t,I)
ylabel('Current, I, A')
%ylim([-4000,8000])

V_t = V_0*(-a/wd).*exp(a.*t).*sin(wd.*t) + V_0*exp(a.*t).*cos(wd.*t);
figure(2)
plot(t,V_t)
ylabel('Voltage, V')

%end