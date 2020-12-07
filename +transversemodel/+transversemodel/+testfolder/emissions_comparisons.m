%% Emissions comparisons version 1

close all
clear all
import transversemodel.subfunctions.*;
% General Parametrs
e = 1.602176634e-19;
m_i = 1.67262192369e-27;

% Test FEE
W = 4.5; %Workfunction 
E_F = 10; % Fermi level
corF = 0.65;  % Develop function for later

% Test TEE
A_G = 1.2016e+06;
T_w= 700;   %Kelvin 

% Test SEE
Te = 1*e;
pow = 22;
n =  10^pow;           % Electron density  [m^-3]
ui0 = sqrt(Te/m_i);      % Ion sheath boundary velocity
E_i = 13.6;             % Ionization energy of hydrogen in electronvolt [eV]
gi = n*ui0; % Ion sheath flux

%% Comparison for Electric field
E_w1 = linspace(0,40)*10^8; % Electric field

[TEE_p1, SEE_p1, FEE_p1] = iter_ems(T_w, W, gi, E_w1, A_G, e, E_i, E_F, corF);

figure(1)
semilogy(E_w1,TEE_p1,E_w1,SEE_p1,E_w1,FEE_p1)
ulm = ceil(log10(max(max([TEE_p1,FEE_p1,SEE_p1]))));   % To set a boundary of interest for the plot
ylim([10^(-ulm) 10^(ulm)])
title(sprintf('Electron emissions (T_W = %u k, n = 10^{%.2e} m^{-3})',T_w,pow))
xlabel('Electric field [V/m]')
ylabel('Current density [A/m^2]')
legend("TEE","SEE","FEE")

%Intersections:

negF = find(FEE_p1>=10^(-ulm),1)         % First non-negligbly small FEE

int_TS = find(min(abs(TEE_p1-SEE_p1))==abs(TEE_p1-SEE_p1))
int_FS = find(min(abs(FEE_p1-SEE_p1))==abs(FEE_p1-SEE_p1))
%int_TF = find((min(abs(TEE_p1-FEE_p1))==abs(TEE_p1-FEE_p1))
int_TF = find(min(abs(TEE_p1(negF:end)-FEE_p1(negF:end)))==abs(TEE_p1-FEE_p1))

iter_E_W = sort([int_TS int_FS int_TF]);

%% Comparison for wall temperature
% Make free plots aroudn the situation detected


T_w2 = linspace(T_w-500,T_w+500);


for j = 1:3
E_w2 = E_w1(iter_E_W(j));
[TEE_p2, SEE_p2, FEE_p2] = iter_ems(T_w2, W, gi, E_w2, A_G, e, E_i, E_F, corF);


figure(1+j)
semilogy(T_w2,TEE_p2,T_w2,SEE_p2,T_w2,FEE_p2)
ulm = ceil(log10(max(max([TEE_p2,FEE_p2,SEE_p2]))));   % To set a boundary of interest for the plot
ylim([10^(-ulm) 10^(ulm)])
title(sprintf('Electron emissions (E_w = %.2e)',E_w2))
xlabel('Wall temperature ')
ylabel('Current density [A/m^2]')
legend("TEE","SEE","FEE")

end
%% Comparison for density
n = linspace(0,4*10^23);

T_w3 = 700; % assumed T_w of 700 K


ui0 =sqrt(Te/m_i);      % Ion sheath boundary velocity
gi = n*ui0; % Ion sheath flux

for j = 1:3
E_w3 = E_w1(iter_E_W(j));
[TEE_p3, SEE_p3, FEE_p3] = iter_ems(T_w3, W, gi, E_w3, A_G, e, E_i, E_F, corF);
figure(4+j)
semilogy(n,TEE_p3,n,SEE_p3,n,FEE_p3)
ulm = ceil(log10(max(max([TEE_p3,FEE_p3,SEE_p3]))));   % To set a boundary of interest for the plot
ylim([10^(-ulm) 10^(ulm)])
title(sprintf('Electron emissions (T_w = 700K, E_w = %.2e)',E_w3))
xlabel('Electron density')
ylabel('Current density [A/m^2]')
legend("TEE","SEE","FEE")
end

%% Iterating function

function [TEE_it, SEE_it, FEE_it] = iter_ems(T_w, W, gi, E_w, A_G, e, E_i, E_F, corF)
import transversemodel.subfunctions.*;
str=[];
if sum([length(T_w) length(E_w) length(gi)]>1)>1
    warning("There are multiple variables introduced as iterables")
    prompt = 'Do you want to continue? Y/N [Y]: ';
    str = input(prompt,'s');
   
end

 
if isempty(str)
   str = 'Y'; 
end

if str=='Y'
    
TEE_it = zeros(100,1);
SEE_it = zeros(100,1);
FEE_it = zeros(100,1);
    
for i= 1:length(T_w)
    for ii = 1:length(E_w)
        for i3 = 1:length(gi)
            TEE_it(i+ii+i3-2) = schottky(T_w(i), W, E_w(ii), A_G)*e;
            SEE_it(i+ii+i3-2) = SEE(gi(i3), E_i, W)*e; %Emissions in A/m^2
            FEE_it(i+ii+i3-2) = FEE(W,E_F,E_w(ii))*e; %Emissions in A/m^2
        end
    end
end

else 
    warning('loop aborted')
end

end