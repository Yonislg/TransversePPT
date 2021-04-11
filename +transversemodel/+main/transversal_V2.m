% Transversal model
% Created 08/04/2020 by Yonis le Grand
clear all
close all
%Physical constants
global k e eps_0 hbar mu_0 m_e imeq;
import transversemodel.subfunctions.*;

k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;
m_i = 1.67262192369e-27;%*17;
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
F_SEE=0; %Include SEE 
if F_SEE
    fprintf('\nSecondary electron emmissions are included');
    SEE_inc = 'yes';
else
    fprintf('\nSecondary electron emmissions are not included');
    SEE_inc = 'no';
end

F_TEE=1; %Include SEE 
if F_TEE
    fprintf('\nThermal electron emmissions are included');
    TEE_inc = 'yes';
else
    fprintf('\nThermal electron emmissions are not included');
    TEE_inc = 'no';
end

F_FEE=0; %Include FEE 
if F_FEE
    fprintf('\nField electron emmissions are included');
    FEE_inc = 'yes';
else
    fprintf('\nField electron emmissions are not included');
    FEE_inc = 'no';
end

%Electrode temepratures
T_wk = 700;           % Electrode temperature (Dummy) [k]
T_wka = 500;          % Anode temperture [k]
T_wkc = 3200;          % Cathode temperture [k]

% Electrode dimension
h = 0.05;               % Distance between electrodes
L = 0.08;               % Length of electrodes

% Material properties
W= 4.4;                 % Workfunction of copper (4.4) in eV according to Frid
A_G1 = 80*10^4;%      % Material Constant 
A_G = 80*10^4*F_TEE%0;%      % Material Constant
E_i = 13.6;             % Ionization energy of hydrogen in electronvolt [eV]
E_F = 7;                % Fermi energy [eV]
%7.7;                   % Ionization energy of copper in electronvolt

% External potential
phi_A =  1000;          % Anode potential 
phi_C = 0;%0;           % Cathode potenital

% Electron properties            
nb = 5*10^22;              % bulk electron density
Te = 2*e;               % Electron Temperature in joules (not eV!)

ne_0 = nb/2;             % Electron density at sheath edge [m^-3]

%ui0 =sqrt(Te/m_i);      % Ion sheath boundary velocity

% Axial inputs
u_ze = 10^4;             % Downstream (axial) flow velocity in m/s
d_J = 0;%.02;               % thickness os the plasma sheet


% neutral density and ion density
Z = 1;                  % ionisation number
a_iz = 1;              % ionisation degree
n_n = (1-a_iz)/a_iz*ne_0;
%ni_b = nb/Z;
ni_0 = ne_0/Z;

%gi = ni_0*sqrt(Te/m_i); % Ion sheath flux



plasma_properties = {Te, ne_0, nb, n_n, Z, m_i};
design_parameters = {T_wka, T_wkc, E_i, A_G, h, L, W, E_F};

%% Initial Guesses for electrode emissions, Wall cleaelectric field
C_guess = 5;

varphi_sf = 0.5*log(2*pi*m_e/m_i) % Guess for sheath potential drop of a floating wall assuming Ti=0

VA_guess= -4.6246*e/Te
%log(2*exp(varphi_sf)/(1+exp(e*(phi_C-phi_A)/Te)));%varphi_sf%%+2
VC_guess=  -14.997*e/Te
%(varphi_sf-C_guess);

%varphi_sf = 0.5*log(2*pi*m_e/m_i) % Guess for sheath potential drop assuming Ti=0


%VA_guess= log(2*exp(varphi_sf)/(1+exp(e*(phi_C-phi_A)/Te)))%varphi_sf%%+2
%VC_guess= phi_C-(phi_A+varphi_sf)%(phi_A-phi_C)-VA_guess%log(2*exp(varphi_sf)/(1+exp(e*(phi_A-phi_C)/Te)))%varphi_sf-10%VA_guess%-4
%gem_guess =  -SEE(gi, E_i, W);% 
%ge_SEE= SEE(gi, E_i, W);  % SEE portion of electron emissions


%% Function solver

% Iterating over different values for wall electric field
    initial_state = init_guessor(VC_guess, VA_guess,T_wka, T_wkc,  Te, nb, ne_0, ni_0, m_i, E_i, E_F, A_G, W, F_SEE, F_FEE)

    %initial_state = [VC_guess,VA_guess,1.76641E+25, 5.64028E-13, 61348447.25, 27478193.15, -7073.859241]
    
    %%
    [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,By,x,fx, jacobian] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C, u_ze, d_J, F_SEE, F_FEE)
    

%%
%import transversemodel.subfunctions.*;
    geC_bolz = ge_bolz(ne_0, Te, -V_C*Te/e)*e
    geC_em*e
    SEE_C = SEE(ni_0*sqrt(Te/m_i), E_i, W)*e
    TEE_C = schottky(T_wkc, W, E_wc, A_G1)*e
    J = -e*nb*uxe
%    FEE_C = FEE(W,E_F, E_wc)
     geA_bolz = ge_bolz(ne_0, Te, -V_A*Te/e)*e
     geA_em*e
     SEE_A = SEE(ni_0*sqrt(Te/m_i), E_i, W)*e
     TEE_A = schottky(T_wka, W, E_wc, A_G1)*e
%     FEE_A = FEE(W,E_F, E_wa)
    
%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx, jacobian] = currentsolver(plasma_properties,design_parameters, x, phi_A, phi_C);
%% Setting up table
InpuValues = {nb, Te/e,phi_A, T_wkc, T_wka, h, d_J, SEE_inc, TEE_inc, FEE_inc};
InputVarNames = {'n_e [m^{-3}]', 'T_e [eV]','Phi_A [V]' , 'T_{w,A} [K]', 'T_{w,c} [K]', 'h [m]','\delta [m]','SEE','TEE','FEE'};
varTypes=[repmat({'double'},1,length(InputVarNames)-3) repmat({'string'},1,3)];
sz= [1 length(InputVarNames)];
InpuTable = table('Size', sz,'VariableTypes',varTypes, 'VariableNames', InputVarNames);
InpuTable(1,:) = InpuValues

OutputVarNames = {'U_{CD} [V]','U_{AB} [V]', 'g*_C [m^2/s]','g*_A [m^2/s]', 'E_{w,c} [V/m]','E_{w,A} [V/m]','u_ze [m/s]','Phi_B', 'Phi_D','B_y [T]'};
varTypes2=repmat({'double'},1,length(OutputVarNames));
sz2= [1 length(OutputVarNames)];
%IniTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);
OutpuTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);
OutpuTable(1,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By}

PotVarNames = {'U_{CD} [V]','U_{AB} [V]','E_{w,c} [V/m]','E_{w,A} [V/m]','Phi_B [V]', 'Phi_D [V]','E_x [V/m]'};
varTypes2 = repmat({'double'},1,length(PotVarNames));
sz2= [1 length(PotVarNames)];
%IniTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);
PoTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', PotVarNames);
PoTable(1,:) = {V_C, V_A, E_wc, E_wa, phi_B, phi_D, (phi_B - phi_D) / h}

EmissioNames = {'g_C [a/m^2]','g*_C [a/m^2]', 'SEE_C [a/m^2]', 'TEE_C [a/m^2]','g_A [a/m^2]','g*_A [a/m^2]', 'SEE_A [a/m^2]', 'TEE_A [a/m^2]', 'J [a/m^2]'};
varTypes3=repmat({'double'},1,length(EmissioNames));
sz3= [1 length(EmissioNames)];
%IniTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);
EmissionTable = table('Size', sz3, 'VariableTypes',varTypes3, 'VariableNames', EmissioNames);
EmissionTable(1,:) = {geC_bolz, geC_em*e, SEE_C, TEE_C, geA_bolz, geA_em*e, SEE_A, TEE_A, J}

ErrorNames = {'E_{w,c} [V/m]','g*_C [m^2/s]','particle conservation','electron continuity', 'g*_A [m^2/s]', 'E_{w,A} [V/m]','momentum balance'};
varTypes2=repmat({'double'},1,length(ErrorNames));
sz2= [1 length(ErrorNames)];
%IniTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);
ErrorTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', ErrorNames);
ErrorTable(1,:) = num2cell(fx)
 

%% functions
function [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,By,x,fx,jacobian] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C, u_ze, d_J, F_SEE, F_FEE)
    global e imeq mu_0;
    import transversemodel.main.total_current;
    japat=spones( [1     0     1     0     1     0     0;
   
     0     0     1     0     1     0     0;
     1     1     1     1     0     0     1;
     1     1     1     1     0     0     0;
     0     0     0     1     0     1     0;
     0     1     0     1     0     1     0;
     1     1     0     0     0     0     1]);
 function stop = outfun(x, optimValues, state)
    stop = false;
    
 end
    options = optimoptions('fsolve','MaxFunctionEvaluations',5e6,'MaxIterations',2e5,'Display','iter','JacobPattern', japat,'PlotFcn',@optimplotfirstorderopt); % @optimplotx,'OutputFcn',@outfun); % 
    Te = plasma_properties{1};
    nb = plasma_properties{3};
    fun =  @(x)total_current(x, plasma_properties,design_parameters, phi_A, phi_C,u_ze, d_J, F_SEE, F_FEE)
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
    By = 0.5*mu_0*e*nb*x(7)*d_J;
end

%% Initial state Guessor
function initial_state = init_guessor(VC_guess, VA_guess,T_wka, T_wkc,  Te, nb, ne_0, ni_0, m_i, E_i, E_F, A_G, W, F_SEE, F_FEE)
global imeq e m_e;
import transversemodel.subfunctions.*;

[EC_guess, gemC_guess] = wall_guess(T_wkc, VC_guess, Te, ne_0, ni_0, m_i, E_i,E_F, A_G, W, F_SEE, F_FEE);
[EA_guess, gemA_guess] = wall_guess(T_wka, VA_guess, Te, ne_0, ni_0, m_i, E_i, E_F, A_G, W, F_SEE, F_FEE);

ue_guess = (gemA_guess -  ge_bolz(ne_0, Te, -VA_guess*Te/e)  - gemC_guess +ge_bolz(ne_0, Te, -VC_guess*Te/e))/2/nb;
%fun = @(x)-e*nb*(1000-VA_guess*Te/e-(0-VC_guess*Te/e))/0.05 - x*m_e*nb*(collRate_ei(nb,1,Te));
%display('init guess')
%[ue_guess,fx] = fsolve(fun,ue_guess1)
%(ge_bolz(ne_0, Te, -VC_guess*Te/e) - gemC_guess)/ne_0-sqrt(Te/m_i)
%sqrt(Te/m_i)+(gemA_guess -  ge_bolz(ne_0, Te, -VA_guess*Te/e))/ne_0
% Collecting initial guesses for variables in one state

if imeq
    initial_state = [VC_guess,VA_guess,gemC_guess,gemA_guess,EC_guess,EA_guess,ue_guess];
else
    initial_state = [VC_guess,VA_guess,gemC_guess,gemA_guess,EC_guess,EA_guess];
end

end

function [E_guess, gem_guess] = wall_guess(T_wk, varphi_sf, Te, ne_0, ni_0, m_i, E_i, E_Fin, A_G, W, F_SEE, F_FEE)
%global;
import transversemodel.subfunctions.*;
if F_SEE
    ge_SEE = SEE(ni_0*sqrt(Te/m_i), E_i, W);
else 
    ge_SEE = 0;
end

if F_FEE
    E_F= E_Fin;
else 
    E_F = 0;
end

E_guess1 = wall_e_field(T_wk, varphi_sf,0,ne_0, Te)
gem_guess1 = schottky(T_wk, W, E_guess1, A_G)+ ge_SEE;

E_guess2 = wall_e_field(T_wk, varphi_sf,gem_guess1,ne_0, Te)
gem_guess2 = schottky(T_wk, W, E_guess2, A_G)+ ge_SEE+ FEE(W,E_F, E_guess2)

E_guess3 = wall_e_field(T_wk, varphi_sf,gem_guess2 ,ne_0, Te)
gem_guess3 = schottky(T_wk, W, E_guess3, A_G)+ ge_SEE+ FEE(W,E_F, E_guess3)

E_guess = wall_e_field(T_wk, varphi_sf,gem_guess3 ,ne_0, Te)
gem_guess = schottky(T_wk, W, E_guess, A_G)+ ge_SEE+ FEE(W,E_F, E_guess)


end