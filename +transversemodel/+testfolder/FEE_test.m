% General Parametrs
import transversemodel.subfunctions.*;
E_W = linspace(0,50)*10^8; 
e = 1.602176634e-19;

% Test FEE
W = 4.5; %Workfunction 
E_F = 10; % Fermi level

%ge = SEE(gi, E_i, W)
for i= 1:100
    em_Fl(i)=FEE(4,7,E_W(i),1); % emissions in m^-3/s
   em_A(i) = log10(em_Fl(i)*e/10^4); %Emissions in A/cm^2
end
plot(E_W,em_A)
%semilogy(E_W,em_A)
title('Verification of FEE function (W = 4.5 eV, E_F=10 eV, T_w \approx 0)')
xlabel('Electric field [V/m]')
ylabel('Common logarithm of current density [A/cm^2]')
ylim([-20 8])