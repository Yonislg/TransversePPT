% Test SEE
% Penning mechanism is valid for energies lower than 1 keV
%test gamma for krypton Ar and tungsten p 514 Frid
import transversemodel.subfunctions.*;
gi  = 1;
E_i = 20;
W = linspace(0,10); % work function of tungsten.
gamma= zeros(100,1);

%ge = SEE(gi, E_i, W)
for i= 1:100
    gamma(i) = SEE(gi, E_i, W(i));
end
plot(W,gamma)
title('Verification of SEE function (E_i=20 eV)')
xlabel('Work function [eV]')
ylabel('\gamma_{se}')
% Also test for copper, find exp data

% Think of how to validate this equation
% is it valid for this range and situation?