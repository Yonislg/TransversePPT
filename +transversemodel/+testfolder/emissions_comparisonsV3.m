%% Version 3 of emission comparison tests


close all
clear all

import transversemodel.subfunctions.*;

% General Parametrs
e = 1.602176634e-19;
pre = [1 12 12 19 19 19]

% Test FEE
W = 4.0; %Workfunction 
E_F = 7; % Fermi level
CorF = 0.7;  % Guess best case scenario

% Test TEE
A_G = 1.2016e+06/2;
Twk = 1000;   %Kelvin 

% Test SEE
E_i = [13.6 11.3 24.4 47.9 17.4 35 62.7] ;                    % Ionization energy of hydrogen in electronvolt [eV]


% Main iterables
iter= 200;
Te = 3*e;
%n =  linspace(10^5,3*10^23,iter);       % Electron density  [m^-3]
n =  logspace(1,25,iter);                % Electron density  [m^-3]
varphi = linspace(-1,-2000,iter);
[N,V] = meshgrid(n,varphi)


%% Comparison for Electric field
%E_w1 = linspace(0,40)*10^8; % Electric field

% GET Inspiration from wall field!
for i = 3:6
m_i = pre(i)*1.67262192369e-27;

[SEE_1, FEE_1, TEE_1, E] = iterfield(Twk, V, N, Te, A_G, E_i(i), W, E_F, CorF,m_i);

SEE_o(i,:) = min(SEE_1);
FEE_o(i,:) = max(TEE_1);
TEE_o(i,:) = max(FEE_1);
E_o(i,:) = max(E);%,[],1);

end
%% Plot
figure(1)
loglog(n,E_o)
title("Density and Electric field")  

figure(2)
hold on
%loglog(n,SEE_o,'-',n,TEE_o,':',n,FEE_o,'--')% or use loglog
plot(n, SEE_o,'-',n,TEE_o,':',n,FEE_o,'--')
%ulm = ceil(log10(max(max([TEE_o,FEE_o,SEE_o]))));   % To set a boundary of interest for the plot
%%ylim([10^(-10) 10^(ulm+5)])
title(sprintf('Maximum Electron emissions (T_e = %u eV)',Te/e))
xlabel('Density [m^{-3}]')
ylabel('Current density [A/m^2]')
%legend("SEE","TEE","FEE")


%% Required functions

function [SEE_o, FEE_o, TEE_o, E] = iterfield(Twk, varphi, n, Te, A_G, E_i, W, E_F, CorF, m_i)
import transversemodel.subfunctions.*;
    % general parameter
    e = 1.602176634e-19;
    
    % Ion velocity
    ui0 = sqrt(Te/m_i);             % Ion sheath boundary velocity
    gi = n.*ui0;                    % Ion sheath flux
    
    % Iterating over the E_w field
    E_w = wall_e_field(Twk, varphi, SEE(gi, E_i, W), n, Te);
    %E = E_w;
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    geT = schottky(Twk, W, E_w, A_G);
    geS = SEE(gi, E_i, W);
    geF = FEE(W,E_F, E_w, CorF);
%     ge(:,1:3);
%     geT(:,1:3);
    %geS(:,1:3)
    %geF(:,1:3)
     %schottky(Twk, W, E_w, A_G) 
     %FEE(W,E_F, E_w, CorF)
     
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
     %schottky(Twk, W, E_w, A_G)
     %FEE(W,E_F, E_w, CorF)
     
    %richtig: 
    E = wall_e_field(Twk, varphi, ge, n, Te);
    
    
    %E(imag(E)~=0) = nan;
    SEE_o = SEE(gi, E_i, W)*e;
    FEE_o = FEE(W,E_F, E_w, CorF)*e;
    TEE_o = schottky(Twk, W, E_w, A_G)*e;
end