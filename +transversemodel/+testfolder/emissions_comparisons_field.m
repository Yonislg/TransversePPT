%% Version 3 of emission comparison tests


close all
clear all

import transversemodel.subfunctions.*;

% General Parametrs
e = 1.602176634e-19;
m_i = 28*1.67262192369e-27;

% Test FEE
W = 4.5; %Workfunction 
E_F = 10; % Fermi energy of copper
CorF = 0.9;  % Guess best case scenario

% Test TEE
A_G = (1.2016e+06);
Twk = 3000;   %Kelvin 

% Test SEE
E_i = 14;                    % Ionization energy of hydrogen in electronvolt [eV]


% Main iterables
%= 200;
Te = 3*e;
E_w = linspace(0,50)*10^8; 


%% Comparison for Electric field


figure(1)
g1 = iterfield(300, E_w, A_G, W, E_F, CorF);
g2 = iterfield(1000, E_w, A_G, W, E_F, CorF);
g3 = iterfield(1700, E_w, A_G, W, E_F, CorF);

semilogy(E_w,g1,'b-',E_w,g2,'r:',E_w,g3,'g--')

%ulm = ceil(log10(max(max([TEE_o,FEE_o,SEE_o]))));   % To set a boundary of interest for the plot
ylim([10^(-20) 10^(8)])
   % xlim([10^(19) 10^(25)])
%title(sprintf('Maximum Electron emissions (T_e = %u eV, T_w = %u K)',Te/e, Twk))
xlabel('Electric field [V/m]')
ylabel('Common logarithm of current density [A/cm^2]')
legend("300 K","1000","1700")


function ge = iterfield(Twk, E_w, A_G, W, E_F, CorF)
import transversemodel.subfunctions.*;
e = 1.602176634e-19;
    ge =  (schottky(Twk, W, E_w, A_G) + FEE(W,E_F, E_w))*e/10^4;
    
end