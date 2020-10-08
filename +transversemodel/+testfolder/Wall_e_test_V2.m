% Wall E Test
% Yonis le Grand
% Tests the influence of density and emission mechanisms on the wall
% electric field.

%% Fixed parameters
close all
clear all
global e m_i
k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;
m_i = 1.67262192369e-27;
mu_0 = 1.2566370614359173e-06;

%% Input parameters
CorF = 0.6;         % Exponentail facors
Twk = 700;          % Wall Electric temperature in Kelvin
Twj  = Twk*k;       % Wall Temperature in Joules
E_i = 13.6;         % Ionization energy of hydrogen in electronvolt [eV]
W = 4.0;            % Workfunction
A_G = 80*10^4;      % TEE material constant
E_F = 10;           % Fermi Energy

%% Iteratble Input parameters
iter= 100;
n = linspace(10,10^24,iter);          % Electron/Ion Density
varphi = linspace(-1,-1000,iter);
% vary over n, varphi, TE
% Compare iterations using with and without FEE, TEE

% make somethingvoer zichtelijk. Just like the older version. Clear
% descriptionofwhen these paramters are relevant. 

% Stars for different Te?
% Use different ways for naming plots so that you can iterate


[N,V] = meshgrid(n,varphi);

C = N;%N.*V;

[E_w, E_w0, E_w1, E_w2] = itertemp(1, N, V, E_i, W, CorF, Twk, A_G, E_F);

fiure(1)
surf(N,V,reldif(E_w,E_w0),C)
title("Relative influence of all emissionsons on E_W")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")
figure(2)
surf(N,V,(E_w-E_w1)./E_w1,C)
title("Relative influence of FEE on E_W compared to SEE+TEE")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")
figure(3)
surf(N,V,(E_w-E_w2)./E_w2,C)
title("Relative influence of SEE only on E_W compared to all")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")

cc=jet(iter);

figure(4)
h1 = plot(n,((E_w-E_w0)./E_w0)');
set(h1, {'color'},num2cell(cc,2))
title("Relative influence of all emissionsons on E_W")
xlabel("density [m^{-3}]")
%ylabel("\Phi[V]")
figure(5)
h2 = plot(n,((E_w-E_w1)./E_w1)');
set(h2, {'color'},num2cell(cc,2))
title("Relative influence of FEE on E_W compared to SEE+TEE")
xlabel("density [m^{-3}]")
%ylabel("\Phi[V]")
figure(6)
h3 = plot(n,((E_w-E_w2)./E_w2)');
xlim([10, 2*10^22]);
title("Relative influence of SEE only on E_W compared to all emissions")
set(h3, {'color'},num2cell(cc,2))
xlabel("density [m^{-3}]")
%ylabel("\Phi[V]")

figure(7)
h4 = plot(varphi,reldif(E_w,E_w0));
title("Relative influence of all emissionsons on E_W")
set(h4, {'color'},num2cell(cc,2))
%xlabel("density [m^{-3}]")
xlabel("\Phi[V]")
figure(8)
h5 = plot(varphi,reldif(E_w,E_w2));

title("Relative influence of SEE only on E_W compared to all emissions")
set(h5, {'color'},num2cell(cc,2))
xlabel("\Phi[V]")
% %%
% varfs =cellfun(@(c) sprintf('%0.1e',c),num2cell([varphi(1) varphi(10:10:100)]),'UniformOutput',false);
% 
% for i = 4:6
%     figure(i)
%     colormap(cc)
%     hc= colorbar;
%     set(hc,'Ticklabels',varfs,'Ticks',linspace(0,1,11))
%     set(hc.Label,'String',"\Phi")%[e*V/T_e]")%'Density in m^{-3}')
%     
% end
% 
% dens =cellfun(@(c) sprintf('%0.1e',c),num2cell([n(1) n(10:10:100)]),'UniformOutput',false);
% 
% for i = 7:8
%     figure(i)
%     colormap(cc)
%     hc= colorbar;
%     set(hc,'Ticklabels',dens,'Ticks',linspace(0,1,11))
%     set(hc.Label,'String','n [m^{-3}]')
%     
% end


%% supporting functions
function diff = reldif(A,B)
    diff1 = (A-B)./B;
    diff1(isnan(A)) = nan;
    diff1(isnan(B)) = nan;
    diff = diff1;
end

function [E_w_comp, E_w_orig, E_w_nF, E_w_nFT] = itertemp(it, n, varphi, E_i, W, CorF, Twk, A_G, E_F)
global e m_i

%for i = 3:it
    Te = it*e;
    
    gi = n.*sqrt(Te/m_i);
    E_w_comp = iterfield(Twk, varphi, gi, n, Te, A_G, E_i, W, E_F, CorF); % with all emission types
    E_w_orig =  iterfield(Twk, varphi, 0, n, Te, 0, E_i, W, 0, CorF); % without any emissions
    E_w_nF = iterfield(Twk, varphi, gi, n, Te, A_G, E_i, W, 0, CorF); % without field emissions
    E_w_nFT = iterfield(Twk, varphi, gi, n, Te, 0, E_i, W, 0, CorF);  % with only SEE
%end

end

function E = iterfield(Twk, varphi, gi, n, Te, A_G, E_i, W, E_F, CorF)

    ge_s = SEE(gi, E_i, W);
    E_w = wall_e_field(Twk, varphi, ge_s, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E = wall_e_field(Twk, varphi, ge, n, Te);
    E(imag(E)~=0) = nan;

end

% [N,V] = meshgrid(n,varphi);
% 
% C = N.*V;
% 
% [E_w, E_w2, ge, ge2] = itertemp(1, N, V, E_i, W, CorF, Twk, A_G, E_F);
% 
% surf(N,V,(E_w-E_w2)./E_w2,C)