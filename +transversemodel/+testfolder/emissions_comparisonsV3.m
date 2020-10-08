%% Version 3 of emission comparison tests


close all
clear all

import transversemodel.subfunctions.*;

% General Parametrs
e = 1.602176634e-19;
m_i = 1.67262192369e-27;

% Test FEE
W = 4.0; %Workfunction 
E_F = 7; % Fermi level
CorF = 0.7;  % Guess best case scenario

% Test TEE
A_G = 1.2016e+06/2;
Twk= 700;   %Kelvin 

% Test SEE
E_i = 13.6;                    % Ionization energy of hydrogen in electronvolt [eV]


% Main iterables
iter= 10;
Te = 3*e;
pow = 23;
n =  linspace(10,3*10^23,iter);       % Electron density  [m^-3]
varphi = linspace(-1,-2000,iter);
[N,V] = meshgrid(n,varphi);


%% Comparison for Electric field
%E_w1 = linspace(0,40)*10^8; % Electric field

% GET Inspiration from wall field!

[SEE_o, FEE_o, TEE_o, E] = iterfield(Twk, V, N, Te, A_G, E_i, W, E_F, CorF);

figure(1)
semilogy(n,E)
title("Density and Electric field")    
figure(2)
semilogy(n, SEE_o','-',n,TEE_o',':',n,FEE_o','--')
ulm = ceil(log10(max(max([TEE_o,FEE_o,SEE_o]))));   % To set a boundary of interest for the plot
ylim([10^(-20) 10^(ulm+5)])
title(sprintf('Electron emissions (T_e = %u eV)',Te/e))
xlabel('Density [V/m]')
ylabel('Current density [A/m^2]')
%legend("SEE","TEE","FEE")



function [SEE_o, FEE_o, TEE_o, E] = iterfield(Twk, varphi, n, Te, A_G, E_i, W, E_F, CorF)
import transversemodel.subfunctions.*;
    % general parameter
    m_i = 1.67262192369e-27;
    e = 1.602176634e-19;
    
    % Ion velocity
    ui0 = sqrt(Te/m_i);             % Ion sheath boundary velocity
    gi = n.*ui0;                    % Ion sheath flux
    
    % Iterating over the E_w field
    E_w = wall_e_field(Twk, varphi, SEE(gi, E_i, W), n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E = wall_e_field(Twk, varphi, ge, n, Te);
    E(imag(E)~=0) = nan;
    SEE_o = SEE(gi, E_i, W)*e;
    FEE_o = FEE(W,E_F, E_w, CorF)*e;
    TEE_o = schottky(Twk, W, E_w, A_G)*e;
end