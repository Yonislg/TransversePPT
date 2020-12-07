% Wall E test
close all
clear all

import transversemodel.subfunctions.*;

k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;
m_i = 1.67262192369e-27;
mu_0 = 1.2566370614359173e-06;
CorF = 0.6;


n= [10^4 10^10 10^24];%linspace(0,10^23);[10^22 10^22 10^22]
%Te = 2*e;
Twk=700;
%ce0 = sqrt(Te./m_e);
Twj  = Twk*k; % Converting to energy units
len = 100;

E_i = 13.6;             % Ionization energy of hydrogen in electronvolt [eV]
W = 4.4;
A_G = 80*10^4;
E_F = 0%10;


varphi = linspace(-0.1,-3);%,len);%(ones(len,1)*linspace(-1,-2000,len))'

%[N,V] = meshgrid(n, varphi)

%[T,V] = meshgrid(T, varphi)

col = 'rgb';
styl = '.:-';

%figure(1)
%hold on
figure(2)
hold on

%%
for i = 1:3
    Te = 3*e;
    gi = n(i).*sqrt(Te/m_i);
    ge_s = SEE(gi, E_i, W);  %ones(len,1)*linspace(10^24,10^27,len)
    [E_w, ge(i,:)] = iterfield(Twk, varphi, ge_s, n(i), Te, A_G, gi, E_i, W, E_F, CorF);
    
    [E_w2, ge2(i,:)] =  iterfield(Twk, varphi, ge_s, n(i), Te, 0, gi, E_i, W, E_F, CorF);%wall_e_field(Twk, varphi, 0, n(i), Te);
    figure(1)
    semilogy(-varphi,E_w,strcat(col(i),'-'),-varphi,E_w2,strcat(col(i),'.'),-varphi,abs(E_w-E_w2),col(i))%,'DisplayName', sprintf('E_w(g_i=SEE), Te = %.u eV', i))
    if i==1
        hold on
    end
    %semilogy(-varphi,E_w2,'.')%,'DisplayName',sprintf('E_w(g_i=0), Te = %.u eV', i))
    %semilogy(-varphi,abs(E_w-E_w2))%,'DisplayName',sprintf('difference, Te = %.u eV', i))
    
    figure(2)
    plot(-varphi,(E_w-E_w2)./E_w2*100,strcat(styl(i),col(i)),'DisplayName',sprintf('n = %.u m^3', n(i)))%,'color',col(i))  'Te = %.u eV', i))%
end

figure(1)
title('Electric Field with and without influence of emissions for')
ylabel('[V/M]')
xlabel('\phi')
%legend

figure(2)
title('Relative difference of Electric Field with and without influence of emissions')
ylabel('%')
xlabel('\phi')
legend

%%
function [E, ge2] = iterfield(Twk, varphi, ge, n, Te, A_G, gi, E_i, W, E_F, CorF)
    import transversemodel.subfunctions.*;
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge2 =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E_w = wall_e_field(Twk, varphi, ge2, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge2 =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E = wall_e_field(Twk, varphi, ge2, n, Te);
    E(imag(E)~=0) = nan;

end

% frac= Te/Twj;
% par= (1-(1-2*varphi*frac).^0.5);
% scnd_term = ge/(n*ce0)*sqrt(Twk/Te).*par;





% [E_w,B,C] = wall_e_field(Twk, V,G,N, Te);
% %A(imag(A)~=0) = nan;
% 
% figure(1)
% surf(N,V,B)
% figure(2)
% surf(N,V,exp(V))
% figure(3)
% surf(N,V,C)
% figure(4)
% surf(N,V,B./C)
% figure(5)
% surf(N,V,B./exp(V))
% figure(6)
% surf(N,V,E_w)

% figure(1)
% plot(ge',A')
% ylabel('Electric Field [V/m]')
% xlabel('Wall Emissions [m^{-2}s^{-1}]')
% colorbar
% figure(2)
% plot(-varphi,A)
% ylabel('Electric Field [V/m]')
% xlabel('Sheath Potential drop [dimensionless (-e*V/T_e)]')
% figure(3)
% semilogx(-scnd_term,A)
% ylabel('Electric Field [V/m]')
% xlabel('Second Term [dimensionless]')
% figure(4)
% surf(ge,varphi,A)
% colorbar
% figure(5)
% contourf(ge,varphi,A)
% colorbar