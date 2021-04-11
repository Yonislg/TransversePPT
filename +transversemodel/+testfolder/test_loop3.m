%% Test Loop 3
% This test loop will iterate over wall temperature vaules, starting at
% the point where the code resolves for the last time

close all
clear all

import transversemodel.subfunctions.*;
import  transversemodel.main.*;
%Physical constants
global k e eps_0 hbar m_e mu_0 imeq;


% To be added :  u_ze 

k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;         
m_i = 1.67262192369e-27;        % Ion mass, assumed to be equal to one proton
mu_0 = 1.2566370614359173e-06;
hbar = 1.0546e-34;

% Flags
imeq=1; %Include momentum balance equation 
if imeq
    fprintf('Momentum equation is included');
else
    fprintf('Momentum equation is not included');
end

%fixed 
W= 4.4; % Workfunction of copper (4.4) in eV
E_i = 13.6;             % Ionization energy of hydrogen in electronvolt [eV]
E_F = 7;                % Fermi energy [eV]
phi_C = 0;           % Cathode potenital
L = 0.08;      
A_G = 80*10^4;          % Material constant for Schottkey equation
h = 0.05;

T_wka = 500;
%T_wkc = 800;

nb = 10^22;
Te = 2*e;
phi_A = 1000
T1 = 960;
T2 = 1100;
r_T_wkc = linspace(T1,T2,101);
dT = (T2 - T1)/(length(r_T_wkc)-1);
u_ze = 10^4;             % Downstream (axial) flow velocity in m/s

d_J = 0;%0.02; i.e. no magnetic field

% Ionisation parameters 
a_iz = 1;     %ionisation degree
Z = 1;          %Ion charge number

iter = 60

C_guess = 5;


%% Setting up table
InputVarNames = {'n_e [m^{-3}]', 'T_e [eV]','Phi_A [V]' , 'T_{w,c} [K]', 'T_{w,A} [K]', 'h [m]','\delta [m]','SEE','TEE','FEE'};
varTypes=[repmat({'double'},1,length(InputVarNames)-3) repmat({'string'},1,3)];
sz= [iter length(InputVarNames)];
Input = table('Size', sz,'VariableTypes',varTypes, 'VariableNames', InputVarNames);

OutputVarNames = {'U_{CD} [V]','U_{AB} [V]', 'g*_C [m^2/s]','g*_A [m^2/s]', 'E_{w,c} [V/m]','E_{w,A} [V/m]','u_ze [m/s]','Phi_B', 'Phi_D','B_y [T]'};
varTypes2=repmat({'double'},1,length(OutputVarNames));
sz2= [iter length(OutputVarNames)];
Output = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);
linit = length(OutputVarNames)-3;
sz4= [iter linit];
varTypes4=repmat({'double'},1,linit);
IniTable = table('Size', sz4, 'VariableTypes',varTypes4, 'VariableNames', OutputVarNames(1:linit));

PotVarNames = {'U_{CD} [V]','U_{AB} [V]','E_{w,c} [V/m]','E_{w,A} [V/m]','Phi_B [V]', 'Phi_D [V]','E_x [V/m]'};
varTypes2 = repmat({'double'},1,length(PotVarNames));
sz2= [iter length(PotVarNames)];
Potential = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', PotVarNames);

EmissioNames = {'g_C [a/m^2]','g*_C [a/m^2]', 'SEE_C [a/m^2]', 'TEE_C [a/m^2]','g_A [a/m^2]','g*_A [a/m^2]', 'SEE_A [a/m^2]', 'TEE_A [a/m^2]', 'J [a/m^2]'};
varTypes3=repmat({'double'},1,length(EmissioNames));
sz3= [iter length(EmissioNames)];
Emissions = table('Size', sz3, 'VariableTypes',varTypes3, 'VariableNames', EmissioNames);

ErrorNames = {'E_{w,c} [V/m]','g*_C [m^2/s]','particle conservation','electron continuity', 'g*_A [m^2/s]', 'E_{w,A} [V/m]','momentum balance','ExitFlag'};
varTypes2=repmat({'double'},1,length(ErrorNames));
sz2= [1 length(ErrorNames)];
Error = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', ErrorNames);



%% Iterations

ctr2 = 0;

F_SEE=0; %Include SEE 
F_TEE=1; %Include TEE
F_FEE=0; %Include FEE 
A_G1 =A_G*F_TEE;%      % Material Constant for shottkey equation

%% First iteration
for te = 1:2
if F_SEE
    SEE_inc = 'yes';
else
    SEE_inc = 'a posteriori';
end

if F_TEE
    TEE_inc = 'yes';
else
    TEE_inc = 'a posteriori';
end

if F_FEE
    FEE_inc = 'yes';
else
    FEE_inc = 'a posteriori';
end
T_wkc = r_T_wkc(te);

plasma_properties = {Te, nb,a_iz,Z, m_i};
design_parameters = {T_wka, T_wkc, E_i, A_G1, h, L, W, E_F};


[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By, x,fx,exitflag,initial_state,output] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C, u_ze, d_J, C_guess,  F_SEE, F_TEE, F_FEE,0);



% results including non resolved
if isreal(x)
    ctr2 = ctr2 + 1;
    InpuValues = {nb, Te/e,phi_A, T_wkc, T_wka, h, d_J, SEE_inc, TEE_inc, FEE_inc};
    Input(ctr2,:) = InpuValues;
    Output(ctr2,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By};
    Potential(ctr2,:) = {V_C, V_A, E_wc, E_wa, phi_B, phi_D, (phi_B - phi_D) / h};
    initial_state = [V_C*e/Te, V_A*e/Te, geC_em, geA_em, E_wc, E_wa, uxe];
    IniTable(ctr2,:) = array2table(initial_state);
    
    geC_bolz = ge_bolz(nb/2, Te, -V_C*Te/e)*e;
    SEE_C = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
    TEE_C = schottky(T_wkc, W, E_wc, A_G)*e;
    J = -e*nb*uxe;
    geA_bolz = ge_bolz(nb/2, Te, -V_A*Te/e)*e;
    SEE_A = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
    TEE_A = schottky(T_wka, W, E_wc, A_G)*e;
    
    Emissions(ctr2,:) = {geC_bolz, geC_em*e, SEE_C, TEE_C, geA_bolz, geA_em*e, SEE_A, TEE_A, J};
    Error(ctr2,:) = array2table([fx exitflag]); %num2cell(fx);
end

end
%% Other iterations


for te = 3:length(r_T_wkc)
    
    T_wkc = r_T_wkc(te);
    %n_n = (1-a_iz)/a_iz*ne_0;
    %By = 0;%.7;%r_B(j);
    
    % Flags
    if F_SEE
        SEE_inc = 'yes';
    else
        SEE_inc = 'a posteriori';
    end
    
    if F_TEE
        TEE_inc = 'yes';
    else
        TEE_inc = 'a posteriori';
    end
    
    if F_FEE
        FEE_inc = 'yes';
    else
        FEE_inc = 'a posteriori';
    end
    
    ne_0=nb/2;

    % neutral density and ion density
    n_n = (1-a_iz)/a_iz*nb;
    %ni_b = nb/Z;
    ni_0 = ne_0/Z;
    
    plasma_properties = {Te, ne_0, nb, n_n, Z, m_i};
    design_parameters = {T_wka, T_wkc, E_i, A_G1, h, L, W, E_F};
   

if F_FEE
    E_F= E_Fin;
else 
    E_F = 0;
end
    
    
    
    [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By ,x,fx, output, exitflag, jacobian] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C,u_ze,d_J,F_SEE, F_FEE);
    
    
    
    % results including non resolved
    if isreal(x)
        ctr2 = ctr2 + 1;
        InpuValues = {nb, Te/e,phi_A, T_wkc, T_wka, h, d_J, SEE_inc, TEE_inc, FEE_inc};
        Input(ctr2,:) = InpuValues;
        
        Output(ctr2,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By};
%         %if i>=4
%             f2 = f1;        
%             f1 = jacobian;
%             initial_state = IniTable{ctr2-1,:} - f1*(IniTable{ctr2-1,:}-IniTable{ctr2-2,:})/(f2-f2);
       %else
            f1 = (IniTable{ctr2-1,:}-IniTable{ctr2-2,:})/dT;
            initial_state = [V_C*e/Te, V_A*e/Te, geC_em, geA_em, E_wc, E_wa, uxe]+f1*dT;
        %end
        %f2 = (IniTable{ctr2-2,:}-IniTable{ctr2-3,:})/dT;%.*IniTable{ctr2-1,:};%-IniTable{ctr2-2,:})./(f1-f2)
        IniTable(ctr2,:) = array2table(initial_state);
        Potential(ctr2,:) = {V_C, V_A, E_wc, E_wa, phi_B, phi_D, (phi_B - phi_D) / h};
        
        
        geC_bolz = ge_bolz(nb/2, Te, -V_C*Te/e)*e;
        SEE_C = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
        TEE_C = schottky(T_wkc, W, E_wc, A_G)*e;
        J = -e*nb*uxe;
        geA_bolz = ge_bolz(nb/2, Te, -V_A*Te/e)*e;
        SEE_A = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
        TEE_A = schottky(T_wka, W, E_wc, A_G)*e;
        
        Emissions(ctr2,:) = {geC_bolz, geC_em*e, SEE_C, TEE_C, geA_bolz, geA_em*e, SEE_A, TEE_A, J};
        Error(ctr2,:) = array2table([fx exitflag]); %num2cell(fx);
    end
    
end


Input(ctr2+1:end,:)=[]
Output(ctr2+1:end,:)=[]
Potential(ctr2+1:end,:)=[]
Emissions(ctr2+1:end,:)=[]
Error(ctr2+1:end,:)=[]
IniTable(ctr2+1:end,:)=[]
%% Plot outputs
for i = 1:7
figure(i)
plot(Input{:,4},Output{:,i})
xlabel(InputVarNames(4))
ylabel(OutputVarNames(i))

end


%% functions

function [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By, x,fx,output,exitflag,jacobian] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C,u_ze,d_J,F_SEE, F_FEE)

    global e imeq mu_0;
    import transversemodel.main.total_current;
    japat=spones( [1     0     1     0     1     0     0;
   
     0     0     1     0     1     0     0;
     1     1     1     1     0     0     1;
     1     1     1     1     0     0     0;
     0     0     0     1     0     1     0;
     0     1     0     1     0     1     0;
     1     1     0     0     0     0     1]);

    options = optimoptions('fsolve','MaxFunctionEvaluations',4.2e4,'MaxIterations',5e3,'Display','none','JacobPattern', japat);%%,'StepTolerance',1e-4);%'PlotFcn',@optimplotfirstorderopt);

    Te = plasma_properties{1};
    nb = plasma_properties{3};
    fun =  @(x)total_current(x, plasma_properties,design_parameters, phi_A, phi_C,u_ze,d_J,F_SEE, F_FEE,0);
    x0 = initial_state; %[varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess,ui0]%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe];
    [x,fx,exitflag,output,jacobian]  = fsolve(fun,x0,options);
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
function initial_state = init_guessor(VC_guess, VA_guess,T_wka, T_wkc,  Te, nb, ne_0, ni_0, m_i, E_i, E_F, A_G, W,F_SEE, F_FEE)
global imeq e ;
import transversemodel.subfunctions.*;

[EC_guess, gemC_guess] = wall_guess(T_wkc, VC_guess, Te, ne_0, ni_0, m_i, E_i,E_F, A_G, W,F_SEE, F_FEE);
[EA_guess, gemA_guess] = wall_guess(T_wka, VA_guess, Te, ne_0, ni_0, m_i, E_i, E_F, A_G, W,F_SEE, F_FEE);

ue_guess = (gemA_guess -  ge_bolz(ne_0, Te, -VA_guess*Te/e)  - gemC_guess +ge_bolz(ne_0, Te, -VC_guess*Te/e))/2/nb;
%(ge_bolz(ne_0, Te, -VC_guess*Te/e) - gemC_guess)/ne_0-sqrt(Te/m_i)
%sqrt(Te/m_i)+(gemA_guess -  ge_bolz(ne_0, Te, -VA_guess*Te/e))/ne_0
% Collecting initial guesses for variables in one state
if imeq
    initial_state = [VC_guess,VA_guess,gemC_guess,gemA_guess,EC_guess,EA_guess,ue_guess];
else
    initial_state = [VC_guess,VA_guess,gemC_guess,gemA_guess,EC_guess,EA_guess];
end

end

function [E_guess, gem_guess] = wall_guess(T_wk, varphi_sf, Te, ne_0, ni_0, m_i, E_i, E_Fin, A_G, W,F_SEE, F_FEE)

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

E_guess1 = wall_e_field(T_wk, varphi_sf,0,ne_0, Te);
gem_guess1 = schottky(T_wk, W, E_guess1, A_G)+ ge_SEE;

E_guess2 = wall_e_field(T_wk, varphi_sf,gem_guess1,ne_0, Te);
gem_guess2 = schottky(T_wk, W, E_guess2, A_G)+ ge_SEE+ FEE(W,E_F, E_guess2);

E_guess3 = wall_e_field(T_wk, varphi_sf,gem_guess2 ,ne_0, Te);
gem_guess3 = schottky(T_wk, W, E_guess3, A_G)+ ge_SEE+ FEE(W,E_F, E_guess3);

E_guess = wall_e_field(T_wk, varphi_sf,gem_guess3 ,ne_0, Te);
gem_guess = schottky(T_wk, W, E_guess, A_G)+ ge_SEE+ FEE(W,E_F, E_guess);


end