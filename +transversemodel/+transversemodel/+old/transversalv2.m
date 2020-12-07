% Transversal model v2
% Testing whether making subfunctions of total_current2 anonymous functions
% improves solver. no apparent difference
% Created 08/04/2020 by Yonis le Grand
%Physical constants
global k e eps_0 hbar m_e m_i mu_0;

k = physconst('boltzmann');
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12;
hbar = 1.0545718176461565e-34;
m_e = 9.1093837015e-31;
m_i = 1.67262192369e-27;
mu_0 = 1.2566370614359173e-06;

W= 4; % Workfunction in eV


% Test values for wall electric field
b = [0 0.8 1.7 2.3 2.8 3.3];    % Values from 
%E_w = b*10^8;

T_wk = 700; % In kelvin (Dummy)
T_wka = 700;
T_wkc = 700;

% T_wk = 500; % In kelvin (Dummy)
% T_w = T_wk*k; % in joules (Dummy)
% T_wka = 500; % Anode temperture in kelvin
% T_wkc = 500; % Cathode temperature in kelvin
% T_wa = T_wka*k; % Anode temperture in joules
% T_wc = T_wkc*k; % Cathode temperature in joules
A_G = 80*10^4;
ne_0 = 10^22;
Te = 1*e; % In joules!! not eV
ui0 =sqrt(Te/m_i);
E_i = 13.6 % Ionization energy of hydrogen in electronvolt
%7.7; % Ionization energy of copper in electronvolt
gi = ne_0*sqrt(Te/m_i);
varphi_sf = 0.5*log(2*pi*m_e/m_i); % Guess for sheath potential drop assuming Ti=0
phi_A =  0; % Anode potential 
phi_C = 1000; % Cathode potenital
h = 0.05;
L = 0.08;

plasma_properties = {Te, ne_0, ui0};
design_parameters = {T_wka, T_wkc, E_i, A_G, h, L, W};


ue_guess = 5000/(e*ne_0*L)
%gem_guess =  -SEE(gi, E_i, W);% 
gi_guess= SEE(ne_0*sqrt(Te/m_i), E_i, W);
gem_guess = schottky(T_wk, W, 0, A_G) %+ SEE(ne_0*sqrt(Te/m_i), E_i, W);

E_guess1 = wall_e_field(T_wk, varphi_sf,gem_guess ,ne_0, Te, ui0);
gem_guess = schottky(T_wk, W, E_guess1, A_G)%+ gi_guess;
E_guess = wall_e_field(T_wk, varphi_sf,gem_guess ,ne_0, Te, ui0)
initial_state = [varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess,ue_guess];

[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C)

function [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D] = currentsolver(plasma_properties,design_parameters, initial_state, phi_A, phi_C)
    global e;
    options = optimoptions('fsolve','MaxFunctionEvaluations',5e4,'MaxIterations',2e4,'Display','iter','PlotFcn',@optimplotfirstorderopt);
    Te = plasma_properties{1};
    fun =  @(x)total_current2(x, plasma_properties,design_parameters, phi_A, phi_C);
    x0 = initial_state %[varphi_sf,varphi_sf,gem_guess,gem_guess,E_guess,E_guess,ui0]%[V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe];
    x = fsolve(fun,x0,options)
    V_C    = x(1);%*Te
    V_A    = x(2);%*Te
    geC_em = x(3);
    geA_em = x(4);
    E_wc   = x(5);
    E_wa   = x(6);
    uxe    = x(7);
    phi_B  = phi_A - V_A*Te/e;
    phi_D  = phi_C + V_C*Te/e;
end
%D=schottky_decrease(E_w)/e
%schottky(T_w, W, E_w, A_G)*10^(-4)

function E_w = wall_e_field(T_wk, varphi_w,g_e,ne_0, Te, ui0)
    global m_i eps_0 m_e k e;
    %varphi_w= varphi_w/Te% already scaled in terms of energy
    T_w = T_wk*k; % Converting to energy units
    ce0 = sqrt(Te/m_e);
    Mi0 = ui0/sqrt(Te/m_i);
    frac= Te/T_w;
    par= (1-sqrt(1-2*varphi_w*frac));
    scnd_term = g_e/(ne_0*ce0)*sqrt(T_w/Te)*par;
    thrd_term = Mi0^2*(sqrt(1-2*varphi_w/Mi0^2) - 1);
    half_dvarphi_dxi_squared = exp(varphi_w)-1 + scnd_term + thrd_term;
    if half_dvarphi_dxi_squared<0
        display(half_dvarphi_dxi_squared)
        display(scnd_term)
        display(thrd_term)
    end
    dvarphi_dxi = sqrt(2*half_dvarphi_dxi_squared);
    E_w = sqrt(ne_0*Te/eps_0)*dvarphi_dxi;
end
function delta_W = schottky_decrease(E_w)
    global e eps_0;
    delta_W = (e^3 * E_w / (4 * pi * eps_0)).^0.5;
end

function ge = schottky(T_w, W, E_w, A_G)
    global e k; 
    W = ones(1,length(E_w))*W*e;
    delta_W = schottky_decrease(E_w);
    %delta_W = [0, 1.07, 1.57, 1.81, 2.01, 2.18]*e
    exponent = -(W - delta_W) /(k *T_w);
    ge = A_G *T_w.^2 / e * exp(exponent)*e; % including the term e depends on parameter A_G
end

% function ge = FEE( W, E_fermi, corr_f, E_w)
%    global hbar m_e e;
%     exponent= -(4*(2*m_e*W).^0.5*W*corr_f)/(3*e*hbar*E_w)
%     ge = e/(4*pi.^2*(W+E_fermi))*(E_fermi/W).^0.5*np.exp(exponent)
% end

function ge = SEE(gi, E_i, W)
    if E_i-W>W
        gamma = 0.016 *(E_i - 2* W);
        ge = gi * gamma;
    else
        ge = 0;
    end
end

function I = ge(ne, Te, V)
    global e m_e;
	I = ne * (Te / (2 * pi * m_e)) .^ 0.5 * exp(-e * V / Te);
end

function jxe_v2 = total_current2(x, plasma_properties,design_parameters, phi_A, phi_C)
    global m_i mu_0 e;
    [Te, ne_0, ui0] = deal(plasma_properties{:});
    [T_wka, T_wkc, E_i, A_G, h, L, W] = deal(design_parameters{:});
    a1= @(a,b)wall_e_field(T_wkc, a,b,ne_0, Te, ui0);
    a2 = @(c)schottky(T_wkc, W, c, A_G);
    a3 = @(d)ge(ne_0, Te, d);
    gi = ne_0*sqrt(Te/m_i);
    jxe_v2(1) = a1(x(1),x(2)) - x(5); % Cathode electric field
    jxe_v2(2) = - x(3) + SEE(gi, E_i, W) + a2(x(5)); %Cathode lectron emissions
    jxe_v2(3) = (x(4) - a3(-x(2)*Te/e)  - x(3) + a3(-x(1)*Te/e))/2 - ne_0*x(7);
    jxe_v2(4) = -2*ne_0*sqrt(Te/m_i) + a3(-x(2)*Te/e)+ a3(-x(1)*Te/e)- x(4)- x(3);
    jxe_v2(5) =  - x(4)  + SEE(gi, E_i, W) + a2(x(6)); % Anode electron emissions
    jxe_v2(6) = a1(x(2),x(4)) -x(6);  %Anode electric field, adapt in future for T_wka =/= T_wkc
    jxe_v2(7) = -(phi_A-x(1)*Te/e-(phi_C+x(2)*Te/e))/h+ 10^3*(e*ne_0*L*x(7)*mu_0); % Bulk current
    
end

% 
% V_C    = x(1)
% V_A    = x(2)
% geC_em = x(3)
% geA_em = x(4)
% E_wc   = x(5)
% E_wa   = x(6)
%uxe    = x(7)

% function jxe_v2 = total_current2(V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe)
%     global m_i W T_w ne_0 Te ui0 A_G;
%     gi = ne_0*sqrt(Te/m_i);
%     jxe_v2(1) = wall_e_field(T_w, V_C,geC_em,ne_0, Te, ui0) -E_wc;
%     jxe_v2(2) = schottky(T_w, W, E_wc, A_G)- geC_em; %+ SEE(gi, E_i, W) 
%     jxe_v2(3) = (geA_em - ge(ne_0, Te, V_A)  - geC_em + ge(ne_0, Te, V_C))/2 - ne_0*uxe;
%     jxe_v2(4) = -2*ne_0*sqrt(Te/m_i) + ge(ne_0, Te, V_A)+ ge(ne_0, Te, V_C)- geA_em- geC_em;
%     jxe_v2(5) = schottky(T_w, W, E_wa, A_G) - geA_em; % + SEE(gi, E_i, W) 
%     jxe_v2(6) = wall_e_field(T_w, V_A,geC_em,ne_0, Te, ui0) -E_wa;    
% end


% 
% function jxe_v1 = total_current(T_w, ne_0, Te, ui0, A_G,gi)
%     global m_i;
%     jxe_v1(1) = wall_e_field(T_w, V_C,geC_em,ne_0, Te, ui0) -E_w;
%     jxe_v1(2) = schottky(T_w, W, E_w, A_G)+ SEE(gi, E_i, W) - geC_em;% -ne_0*sqrt(Te/m_i);
%     jxe_v1(3) = ge(ne_0, Te, V_C) - geC_em -ne_0*sqrt(Te/m_i) - ne_0*uxe;
%     
% 
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