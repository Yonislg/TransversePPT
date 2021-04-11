%% Version 3 of emission comparison tests


%close all
clear all

import transversemodel.subfunctions.*;

% General Parametrs
e = 1.602176634e-19;
m_i = 28*1.67262192369e-27;

% Test FEE
W = 4.5; %Workfunction 
E_F = 7; % Fermi energy of copper
CorF = 0.9;  % Guess best case scenario

% Test TEE
A_G = (1.2016e+06)/2;
Twk = linspace(2500,4000, 20);   %Kelvin 

% Test SEE
E_i = 14;                    % Ionization energy of hydrogen in electronvolt [eV]


% Main iterables
%= 200;
Te = 3*e;
%n =  linspace(10^5,3*10^23,iter);       % Electron density  [m^-3]
n =  logspace(20,24,200);%iter);       % Electron density  [m^-3]
varphi = -logspace(1,3,5);%iter);
[N,T] = meshgrid(n,Twk);


%% Comparison for Electric field
%E_w1 = linspace(0,40)*10^8; % Electric field

% GET Inspiration from wall field!
for i = 1:5
[SEE_o, FEE_o, TEE_o, E] = iterfield(T, varphi(i), N, Te, A_G, E_i, W, E_F, CorF);
cc=jet(20);
tiLabls1 = cellfun(@(c) sprintf('%0.1e',c),num2cell(linspace(2500,4000,11)),'UniformOutput',false)
figure(2*i-1)

hbv = loglog(n,E');
     set(hbv, {'color'},num2cell(cc,2))
         colormap(cc)
     hc= colorbar;
     set(hc,'Ticklabels',tiLabls1,'limit',[0 1],'Ticks',linspace(0,1,11))
%     set(hc.Label,'String','Potential drop in V')
%     xlabel('Density [m^{-3}]')
%     xlim([10^(19) 10^(25)])
ylabel('Electric Field [V/m]')
title(sprintf('Electric field (T_e = %u eV, U = %i )',Te/e, varphi(i)))  
figure(2*i)
hbr = loglog(n,TEE_o,':')
     set(hbr, {'color'},num2cell(cc,2))
         colormap(cc)
     hr= colorbar;
     set(hr,'Ticklabels',tiLabls1,'limit',[0 1],'Ticks',linspace(0,1,11))
     hold on
     loglog(n,SEE_o,'-',n,FEE_o,'--')
%loglog(n, max(SEE_o),'-',n,max(TEE_o),':',n,max(FEE_o),'--')
ulm = ceil(log10(max(max([TEE_o,FEE_o,SEE_o]))));   % To set a boundary of interest for the plot
ylim([10^(-10) 10^(ulm+5)])
   % xlim([10^(19) 10^(25)])
title(sprintf('Maximum Electron emissions (T_e = %u eV, U = %i)',Te/e, varphi(i)))
xlabel('Density [m^{-3}]')
ylabel('Current density [A/m^2]')
legend("SEE","TEE","FEE")

end
function [SEE_o, FEE_o, TEE_o, E] = iterfield(Twk, varphi, n, Te, A_G, E_i, W, E_F, CorF)
import transversemodel.subfunctions.*;
    % general parameter
    m_i = 1.67262192369e-27;
    e = 1.602176634e-19;
    
    % Ion velocity
    ui0 = sqrt(Te/m_i);             % Ion sheath boundary velocity
    gi = n.*ui0/2;                    % Ion sheath flux
    
    E_w = wall_e_field(Twk, varphi, 0, n, Te);
    
    ge =  schottky(Twk, W, E_w, A_G);%+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    
    % Iterating over the E_w field
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W); %+ FEE(W,E_F, E_w, CorF);
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w);
    E = wall_e_field(Twk, varphi, ge, n, Te);
    %E(imag(E)~=0) = nan;
    SEE_o = SEE(gi, E_i, W)*e;
    FEE_o = FEE(W,E_F, E_w, CorF)*e;
    TEE_o = schottky(Twk, W, E_w, A_G)*e;
end