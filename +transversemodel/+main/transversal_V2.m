% Transversal model
% Created 08/04/2020 by Yonis le Grand
clear all
close all
%Physical constants
global k e eps_0 hbar mu_0 m_e u_ze imeq F_SEE F_FEE F_TEE;
import transversemodel.subfunctions.*;

k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;
m_i = 1.67262192369e-27*17;
mu_0 = 1.2566370614359173e-06;
hbar = 1.0546e-34; % J*s

% Flags
imeq=1; %Include momentum balance equation 
if imeq
    fprintf('Momentum equation is included');
else
    fprintf('Momentum equation is not included');
end

% Flags
F_SEE=1; %Include SEE 
if F_SEE
    fprintf('\nSecondary electron emmissions are included');
else
    fprintf('\nSecondary electron emmissions are not included');
end

F_TEE=1; %Include SEE 
if F_TEE
    fprintf('\nThermal electron emmissions are included');
else
    fprintf('\nThermal electron emmissions are not included');
end

F_FEE=0; %Include FEE 
if F_FEE
    fprintf('\nField electron emmissions are included');
else
    fprintf('\nField electron emmissions are not included');
end
        
%u_ze = 10000;

W= 4.4; % Workfunction of copper (4.4) in eV according to Frid

%Electrode temepratures
T_wk = 700;           % Electrode temperature (Dummy) [k]
T_wka = 900;          % Anode temperture [k]
T_wkc = 600;          % Cathode temperture [k]

A_G = 80*10^4*F_TEE%0;%      % Material Constant
%ne_0 =            
nb = 2.7*10^22;              % bulk electron density
ne_0 =nb/2;             % Electron density at sheath edge [m^-3]
Te = 3*e;               % Electron Temperature in joules (not eV!)
ui0 =sqrt(Te/m_i);      % Ion sheath boundary velocity
E_i = 13.6;             % Ionization energy of hydrogen in electronvolt [eV]
E_F = 7;                % Fermi energy [eV]
%7.7;                   % Ionization energy of copper in electronvolt

phi_A =  1200;          % Anode potential 
phi_C = 0;%0;           % Cathode potenital
h = 0.01;               % Distance between electrodes
L = 0.08;               % Length of electrodes
Z = 2;                  % ionisation number
a_iz = 1;              % ionisation degree

% neutral density and ion density
n_n = (1-a_iz)/a_iz*ne_0;
n_i = ne_0/Z;
gi = n_i*sqrt(Te/m_i); % Ion sheath flux


plasma_properties = {Te, ne_0, n_n, Z, ui0};
design_parameters = {T_wka, T_wkc, E_i, A_G, h, L, W, E_F};

%% Initial Guesses for electrode emissions, Wall cleaelectric field
C_guess = 5;

varphi_sf = 0.5*log(2*pi*m_e/m_i); % Guess for sheath potential drop of a floating wall assuming Ti=0

VA_guess= log(2*exp(varphi_sf)/(1+exp(e*(phi_C-phi_A)/Te)));%varphi_sf%%+2
VC_guess=  varphi_sf-C_guess;

%varphi_sf = 0.5*log(2*pi*m_e/m_i) % Guess for sheath potential drop assuming Ti=0


%VA_guess= log(2*exp(varphi_sf)/(1+exp(e*(phi_C-phi_A)/Te)))%varphi_sf%%+2
%VC_guess= phi_C-(phi_A+varphi_sf)%(phi_A-phi_C)-VA_guess%log(2*exp(varphi_sf)/(1+exp(e*(phi_A-phi_C)/Te)))%varphi_sf-10%VA_guess%-4
%gem_guess =  -SEE(gi, E_i, W);% 
%ge_SEE= SEE(gi, E_i, W);  % SEE portion of electron emissions

%% Function solver

% Iterating over different values for wall electric field
    initial_state = init_guessor(VC_guess, VA_guess,T_wka, T_wkc,  Te, ne_0, m_i, E_i, E_F, A_G, W);

    [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx, jacobian] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C)

%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx, jacobian] = currentsolver(plasma_properties,design_parameters, x, phi_A, phi_C);


%% functions
function [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx,jacobian] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C)
    global e imeq;
    import transversemodel.main.total_current;
    japat=spones( [1     0     1     0     1     0     0;
   
     0     0     1     0     1     0     0;
     1     1     1     1     0     0     1;
     1     1     1     1     0     0     0;
     0     0     0     1     0     1     0;
     0     1     0     1     0     1     0;
     1     1     0     0     0     0     1]);
    options = optimoptions('fsolve','MaxFunctionEvaluations',5e4,'MaxIterations',2e3,'Display','iter','JacobPattern', japat);%,'PlotFcn',@optimplotfirstorderopt);
    Te = plasma_properties{1};
    fun =  @(x)total_current(x, plasma_properties,design_parameters, phi_A, phi_C);
    x0 = initial_state %[varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess,ui0]%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe];
    [x,fx,exitflag,output,jacobian]  = fsolve(fun,x0,options)
    V_C    = x(1)*Te/e;%*Te
    V_A    = x(2)*Te/e;%*Te
    geC_em = x(3);
    geA_em = x(4);
    E_wc   = x(5);
    E_wa   = x(6);
    if imeq
        uxe    = x(7);
    else
        uxe    = 0;
    end
    phi_B  = phi_A - V_A;
    phi_D  = phi_C - V_C;
end

%% Initial state Guessor
function initial_state = init_guessor(VC_guess, VA_guess,T_wka, T_wkc,  Te, ne_0, m_i, E_i, E_F, A_G, W)
global imeq e ;
import transversemodel.subfunctions.*;

[EC_guess, gemC_guess] = wall_guess(T_wkc, VC_guess, Te, ne_0, m_i, E_i,E_F, A_G, W);
[EA_guess, gemA_guess] = wall_guess(T_wka, VA_guess, Te, ne_0, m_i, E_i, E_F, A_G, W);

ue_guess = (gemA_guess -  ge_bolz(ne_0, Te, -VA_guess*Te/e)  - gemC_guess +ge_bolz(ne_0, Te, -VC_guess*Te/e))/2/ne_0;
%(ge_bolz(ne_0, Te, -VC_guess*Te/e) - gemC_guess)/ne_0-sqrt(Te/m_i)
%sqrt(Te/m_i)+(gemA_guess -  ge_bolz(ne_0, Te, -VA_guess*Te/e))/ne_0
% Collecting initial guesses for variables in one state
if imeq
    initial_state = [VC_guess,VA_guess,gemC_guess,gemA_guess,EC_guess,EA_guess,ue_guess];
else
    initial_state = [VC_guess,VA_guess,gemC_guess,gemA_guess,EC_guess,EA_guess];
end

end

function [E_guess, gem_guess] = wall_guess(T_wk, varphi_sf, Te, ne_0, m_i, E_i, E_Fin, A_G, W)
global F_SEE F_FEE;
import transversemodel.subfunctions.*;
if F_SEE
    ge_SEE = SEE(ne_0*sqrt(Te/m_i), E_i, W);
else 
    ge_SEE = 0;
end

if F_FEE
    E_F= E_Fin;
else 
    E_F = 0;
end

E_guess1 = wall_e_field(T_wk, varphi_sf,0,ne_0, Te);
gem_guess1 = schottky(T_wk, W, E_guess1, A_G)+ ge_SEE;

E_guess2 = wall_e_field(T_wk, varphi_sf,gem_guess1,ne_0, Te);
gem_guess2 = schottky(T_wk, W, E_guess2, A_G)+ ge_SEE+ FEE(W,E_F, E_guess2);

E_guess3 = wall_e_field(T_wk, varphi_sf,gem_guess2 ,ne_0, Te);
gem_guess3 = schottky(T_wk, W, E_guess3, A_G)+ ge_SEE+ FEE(W,E_F, E_guess3);

E_guess = wall_e_field(T_wk, varphi_sf,gem_guess3 ,ne_0, Te);
gem_guess = schottky(T_wk, W, E_guess, A_G)+ ge_SEE+ FEE(W,E_F, E_guess);


end