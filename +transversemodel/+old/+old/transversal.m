% Transversal model
% Created 08/04/2020 by Yonis le Grand

%Physical constants
global k e eps_0 hbar m_e m_i mu_0 u_ze gi imeq;

k = physconst('boltzmann');
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
hbar = 1.0545718176461565e-34;
m_e = 9.1093837015e-31;
m_i = 1.67262192369e-27;
mu_0 = 1.2566370614359173e-06;

% Flags
imeq=0; %Include momentum equation
if imeq
    fprintf('Momentum equation is included');
else
    fprintf('Momentum equation is not included');
end
        
u_ze = 10000;

W= 4; % Workfunction in eV

T_wk = 700; % In kelvin (Dummy)
T_wka = 700;
T_wkc = 700;

% T_wk = 500; % In kelvin (Dummy)
% T_w = T_wk*k; % in joles (Dummy)
% T_wka = 500; % Anode temperture in kelvin
% T_wkc = 500; % Cathode temperature in kelvin
% T_wa = T_wka*k; % Anode temperture in joules
% T_wc = T_wkc*k; % Cathode temperature in joules
A_G = 0;%80*10^4;
ne_0 = 10^22;
Te = 1*e; % In joules!! not eV
ui0 =sqrt(Te/m_i);
E_i = 13.6 % Ionization energy of hydrogen in electronvolt
%7.7; % Ionization energy of copper in electronvolt
gi = ne_0*sqrt(Te/m_i);
varphi_sf = 0.3*log(2*pi*m_e/m_i); % Guess for sheath potential drop assuming Ti=0
phi_A =  0; % Anode potential 
phi_C = 1000; % Cathode potenital
h = 0.05;
L = 0.08;

gi = ne_0*sqrt(Te/m_i);

plasma_properties = {Te, ne_0, ui0};
design_parameters = {T_wka, T_wkc, E_i, A_G, h, L, W};

%gem_guess =  -SEE(gi, E_i, W);% 
gi_guess= SEE(ne_0*sqrt(Te/m_i), E_i, W);
%gem_guess = schottky(T_wk, W, 0, A_G) %+ SEE(ne_0*sqrt(Te/m_i), E_i, W);
%ue_guess = 1.4*10^26/(ne_0);

E_guess1 = wall_e_field(T_wk, varphi_sf,0 ,ne_0, Te, ui0);
gem_guess = schottky(T_wk, W, E_guess1, A_G)%+ gi_guess;
E_guess = wall_e_field(T_wk, varphi_sf,gem_guess ,ne_0, Te, ui0)
ue_guess = sqrt(Te/m_i)+(-gem_guess +  ge_bolz(ne_0, Te, -varphi_sf*Te/e))/ne_0
%%
if imeq
    initial_state = [varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess,ue_guess];
else
    initial_state = [varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess];
end

[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C);

function [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C)
    global e imeq;
    options = optimoptions('fsolve','MaxFunctionEvaluations',5e4,'MaxIterations',2e3,'Display','iter');%,'PlotFcn',@optimplotfirstorderopt);
    Te = plasma_properties{1};
    fun =  @(x)total_current(x, plasma_properties,design_parameters, phi_A, phi_C);
    x0 = initial_state %[varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess,ui0]%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe];
    [x,fx] = fsolve(fun,x0,options)
    V_C    = x(1);%*Te
    V_A    = x(2);%*Te
    geC_em = x(3);
    geA_em = x(4);
    E_wc   = x(5);
    E_wa   = x(6);
    if imeq
        uxe    = x(7);
    else
        uxe    = 0;
    end
    phi_B  = phi_A - V_A*Te/e;
    phi_D  = phi_C + V_C*Te/e;
end
%D=schottky_decrease(E_w)/e
%schottky(T_w, W, E_w, A_G)*10^(-4)


% function ge = FEE( W, E_fermi, corr_f, E_w)
%    global hbar m_e e;
%     exponent= -(4*(2*m_e*W).^0.5*W*corr_f)/(3*e*hbar*E_w)
%     ge = e/(4*pi.^2*(W+E_fermi))*(E_fermi/W).^0.5*np.exp(exponent)
% end









% 
% function jxe = compute_jxe(z , properties, varphi_A, varphi_C)
%     % z is a vector with the axial positions of the current solution point
%     %properties is an array with the density, axial velocity, electron temperature at the points z
%     ne = properties(1)
%     ux = properties(2)
%     Te = properties(3)
% 
%     %V_C = V_A
%     geA = sheath_current(ne,Te,V_A)
%     geC = sheath_current(ne, Te, V_C)
% 
%     jxe = (geA_em - geA - geC_em+ geC)/2
% end