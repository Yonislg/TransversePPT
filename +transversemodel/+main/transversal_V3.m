% Transversal model V3
% This code prepares inputs for the total_current.m by making an initial
% estimates, it than processes the output.
% This version (V3) can be called upon by other scripts that determine
% inputs and outputs, which is preferable when making a large amount of iterations over wide ranges.
% The other version (V2) is integrated in a testing script which is
% benefificial when making quick simple calculations.
% Created 08/04/2020 by Yonis le Grand



function [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx, exitflag,initial_state,output] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C,u_ze, By, C_guess)%, F_TEE, F_FEE, F_SEE, imeq)


%Physical constants
global k e eps_0 hbar m_e mu_0 F_FEE F_TEE;
import transversemodel.subfunctions.*;

% To be included:  u_ze 

%% Physics constants
k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;
%m_i = 1.67262192369e-27;               
mu_0 = 1.2566370614359173e-06;
hbar = 1.0546e-34;


   
    [Te, nb, a_iz, Z, m_i] = deal(plasma_properties{:});
    ne_0=nb;%/2!!
    [T_wka, T_wkc, E_i, A_G, h, L, W, E_Fin] = deal(design_parameters{:});

    % neutral density and ion density
    n_n = (1-a_iz)/a_iz*nb;
    %ni_b = nb/Z;
    ni_0 = ne_0/Z;
    
    plasma_properties = {Te, ne_0, nb, n_n, Z, m_i};
   

if F_FEE
    E_F= E_Fin;
else 
    E_F = 0;
end
    
    
A_G =A_G*F_TEE;%      % Material Constant for shottkey equation

%% Initial Guesses for Wall Electric field
varphi_sf = 0.5*log(2*pi*m_e/m_i); % Guess for sheath potential drop of a floating wall assuming Ti=0

VA_guess= log(2*exp(varphi_sf)/(1+exp(e*(phi_C-phi_A)/Te)));%varphi_sf%%+2
VC_guess=  varphi_sf-C_guess;% log(2*exp(varphi_sf)/(1+exp(e*(phi_A-phi_C)/Te)))%varphi_sf-10%VA_guess%-4% ( phi_C-(phi_A+varphi_sf))/2;%
%gem_guess =  -SEE(gi, E_i, W);% 
%ge_SEE= SEE(gi, E_i, W);  % SEE portion of electron emissions

%% Function solver

% Iterating over different values for wall electric field
    initial_state = init_guessor(VC_guess, VA_guess,T_wka, T_wkc,  Te, ne_0, ni_0, m_i, E_i, E_F, A_G, W);


    [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx, output, exitflag] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C,u_ze,By);


%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx, jacobian] = currentsolver(plasma_properties,design_parameters, x, phi_A, phi_C);
end

%% functions

function [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx,output,exitflag] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C,u_ze,By)

    global e imeq;
    import transversemodel.main.total_current;
    japat=spones( [1     0     1     0     1     0     0;
   
     0     0     1     0     1     0     0;
     1     1     1     1     0     0     1;
     1     1     1     1     0     0     0;
     0     0     0     1     0     1     0;
     0     1     0     1     0     1     0;
     1     1     0     0     0     0     1]);

    options = optimoptions('fsolve','MaxFunctionEvaluations',4.2e3,'MaxIterations',5e2,'Display','none','JacobPattern', japat,'StepTolerance',1e-4);%'PlotFcn',@optimplotfirstorderopt);

    Te = plasma_properties{1};
    fun =  @(x)total_current(x, plasma_properties,design_parameters, phi_A, phi_C,u_ze,By);
    x0 = initial_state; %[varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess,ui0]%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe];
    [x,fx,exitflag,output]  = fsolve(fun,x0,options);
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
function initial_state = init_guessor(VC_guess, VA_guess,T_wka, T_wkc,  Te, ne_0, ni_0, m_i, E_i, E_F, A_G, W)
global imeq e ;
import transversemodel.subfunctions.*;

[EC_guess, gemC_guess] = wall_guess(T_wkc, VC_guess, Te, ne_0, ni_0, m_i, E_i,E_F, A_G, W);
[EA_guess, gemA_guess] = wall_guess(T_wka, VA_guess, Te, ne_0, ni_0, m_i, E_i, E_F, A_G, W);

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

function [E_guess, gem_guess] = wall_guess(T_wk, varphi_sf, Te, ne_0, ni_0, m_i, E_i, E_Fin, A_G, W)
global F_SEE F_FEE;
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