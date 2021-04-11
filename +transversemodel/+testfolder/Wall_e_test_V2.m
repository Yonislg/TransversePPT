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
CorF = 0.8;         % Exponentail facors
Twk = 3500;          % Wall Electric temperature in Kelvin
Twj  = Twk*k;       % Wall Temperature in Joules
E_i = 13.6;         % Ionization energy of hydrogen in electronvolt [eV]
W = 4.4;            % Workfunction
A_G = 80*10^4;      % TEE material constant
E_F = 10;           % Fermi Energy

%% Iteratble Input parameters
iter= 100;
n = linspace(10^21,10^22,10);          % Electron/Ion Density
varphi = linspace(0.5,400,iter);
% vary over n, varphi, TE
% Compare iterations using with and without FEE, TEE

% make somethingvoer zichtelijk. Just like the older version. Clear
% descriptionofwhen these paramters are relevant. 

% Stars for different Te?
% Use different ways for naming plots so that you can iterate


[N,V] = meshgrid(n,varphi);

C = N;%N.*V;

%[E_w, E_w0, E_w1, E_w2, E_w3] = itertemp(3, N, V, E_i, W, CorF, Twk, A_G, E_F);

[E_all, E_orig, E_nF, E_nFT, E_nSF] = itertemp(2, N, -V, E_i, W, CorF, Twk, A_G, E_F);
%%
figure(1)
surf(N,V,(E_all-E_orig)./E_orig,C)%(N,V,reldif(E_w,E_w0),C)
title("Relative influence of all emissionsons on E_W")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")
figure(2)
surf(N,V,(E_all-E_orig)./E_orig,C)
title("Relative influence of TEE on E_W compared to none")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")
figure(3)
surf(N,V,(E_nFT-E_all)./E_nFT,C)
title("Relative influence of SEE only on E_W compared to all")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")
figure(11)
surf(N,V,E_all,C)%(N,V,reldif(E_w,E_w0),C)
title("E_W including all emissions")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")
figure(12)
surf(N,V,E_orig,C)%(N,V,reldif(E_w,E_w0),C)
title("E_W no emissions")
xlabel("density [m^{-3}]")
ylabel("\Phi[V]")

cc=jet(iter);
%% Plotting for n
figure(4)
h1 = plot(n,((E_all-E_orig)./E_orig)');
set(h1, {'color'},num2cell(cc,2))
title("Relative influence of all emissionsons on E_W")
xlabel("density [m^{-3}]")
%ylabel("\Phi[V]")
figure(5)
h2 = plot(n,((E_all-E_nSF)./E_nSF)');
set(h2, {'color'},num2cell(cc,2))
title("Relative influence of TEE on E_W compared all")
xlabel("density [m^{-3}]")
%ylabel("\Phi[V]")
figure(6)
h3 = plot(n,((E_orig-E_nFT)./E_orig)');
%xlim([10, 2*10^22]);
title("Relative influence of SEE only on E_W compared to none")
set(h3, {'color'},num2cell(cc,2))
xlabel("density [m^{-3}]")
%ylabel("\Phi[V]")

%% Plotting for varphi
figure(7)
h4 = plot(varphi,reldif(E_all,E_orig));
title("Relative influence of all emissionsons on E_W")
set(h4, {'color'},num2cell(cc,2))
%xlabel("density [m^{-3}]")
xlabel("\Phi[V]")

figure(8)
h5 = plot(varphi,reldif(E_all,E_nFT));
title("Relative influence of SEE only on E_W compared to all emissions")
set(h5, {'color'},num2cell(cc,2))
xlabel("\Phi[V]")

figure(9)
h6 = plot(varphi,reldif(E_all,E_nSF));
title("Relative influence of TEE only on E_W compared to all emissions")
set(h6, {'color'},num2cell(cc,2))
xlabel("\Phi[V]")

figure(10)
h7 = plot(varphi,reldif(E_nSF,E_orig));
title("Relative influence of TEE only on E_W compared to none")
set(h7, {'color'},num2cell(cc,2))
xlabel("\Phi[V]")

figure(13)
h8 = plot(varphi,E_all,varphi,E_orig);
title("E_w")
legend("without emissions", "with emissions")
%set(h7, {'color'},num2cell(cc,2))
xlabel("\Phi[V]")
%%
varfs =cellfun(@(c) sprintf('%0.1e',c),num2cell([varphi(1) varphi(10:10:100)]),'UniformOutput',false);

for i = 4:6
    figure(i)
    colormap(cc)
    hc= colorbar;
    set(hc,'Ticklabels',varfs,'Ticks',linspace(0,1,11))
    set(hc.Label,'String',"\Phi")%[e*V/T_e]")%'Density in m^{-3}')
    
end

dens =cellfun(@(c) sprintf('%0.1e',c),num2cell([n(1) n(10:10:100)]),'UniformOutput',false);

for i = 7:10
    figure(i)
    colormap(cc)
    hc= colorbar;
    set(hc,'Ticklabels',dens,'Ticks',linspace(0,1,11))
    set(hc.Label,'String','n [m^{-3}]')
    
end


%% supporting functions
function diff = reldif(A,B)
    diff1 = (A-B)./B;
    diff1(isnan(A)) = nan;
    diff1(isnan(B)) = nan;
    diff = diff1;
end

function [E_all, E_orig, E_nF, E_nFT, E_nSF] = itertemp(it, n, varphi, E_i, W, CorF, Twk, A_G, E_F)
global e m_i

%for i = 3:it
    Te = it*e;
    
    gi = n.*sqrt(Te/m_i);
    E_all = iterfield(Twk, varphi, gi, n, Te, A_G, E_i, W, E_F, CorF); % with all emission types
    E_orig =  iterfield(Twk, varphi, 0, n, Te, 0, E_i, W, 0, CorF); % without any emissions
    E_nF = iterfield(Twk, varphi, gi, n, Te, A_G, E_i, W, 0, CorF); % without field emissions
    E_nFT = iterfield(Twk, varphi, gi, n, Te, 0, E_i, W, 0, CorF);  % with only SEE
    E_nSF = iterfield(Twk, varphi, 0, n, Te, A_G, E_i, W, 0, CorF);  % with only TEE
%end

end

function E = iterfield(Twk, varphi, gi, n, Te, A_G, E_i, W, E_F, CorF)

    import transversemodel.subfunctions.*;

    ge_s = SEE(gi, E_i, W);
    E_w = wall_e_field(Twk, varphi, ge_s, n, Te);
    %E_w(imag(E_w)~=0) = nan;
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    %E_w(imag(E_w)~=0) = nan;
%     E = zeros(length(n));
%     for i = 1:length(n)
%         for j = 1:length(n)
%                 if length(gi)==1
%                     fun = @(x)sysfunc(x,Twk,W,A_G,gi,E_i,E_F,CorF,varphi(i,j),n(i,j),Te);
%                 else
%                     fun = @(x)sysfunc(x,Twk,W,A_G,gi(i,j),E_i,E_F,CorF,varphi(i,j),n(i,j),Te);
%                 end
%                 x0 = [ge(i,j),E_w(i,j)];
%                 x = fsolve(fun,x0);
%                 E(i,j)= x(2);
%         end
%     end
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E = wall_e_field(Twk, varphi, ge, n, Te);
    
    s=0;
    while mean(mean(abs(E-E_w)))<0.001&&s<100
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E_w = wall_e_field(Twk, varphi, ge, n, Te);
    
    ge =  schottky(Twk, W, E_w, A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, E_w, CorF);
    E = wall_e_field(Twk, varphi, ge, n, Te);
    s=s+1;
    end 
    
    E - E_w
    E(imag(E)~=0) = nan;

end

function F1 = sysfunc(x,Twk,W,A_G,gi,E_i,E_F,CorF,varphi,n,Te)
     import transversemodel.subfunctions.*;

        F1(1) = schottky(Twk, W, x(2), A_G)+ SEE(gi, E_i, W) + FEE(W,E_F, x(2), CorF) -x(1);
        F1(2) = wall_e_field(Twk, varphi, x(1), n, Te) - x(2);

end

% [N,V] = meshgrid(n,varphi);
% 
% C = N.*V;
% 
% [E_w, E_w2, ge, ge2] = itertemp(1, N, V, E_i, W, CorF, Twk, A_G, E_F);
% 
% surf(N,V,(E_w-E_w2)./E_w2,C)