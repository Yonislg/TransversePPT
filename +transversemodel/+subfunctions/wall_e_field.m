
function E_w= wall_e_field(T_wk, varphi_w,g_e,ne_0, Te)%, ui0)%[E_w,scnd_term,thrd_term] 
    % Physics constants
    k = physconst('boltzmann');     % Bolzmann Constant [J/K]
    e = 1.602176634e-19;
    eps_0 = 8.854187817620389e-12; 
    m_e = 9.1093837015e-31;
    m_i = 1.67262192369e-27;
    mu_0 = 1.2566370614359173e-06;
    hbar = 1.0546e-34; % J*s
    
    
    %varphi_w= varphi_w/Te% already scaled in terms of energy
    T_w = T_wk*k; % Converting to energy units
    ce0 = sqrt(Te/m_e);
    Mi0 = 1;%ui0/sqrt(Te/m_i);
    frac= Te/T_w;
    par= (1-(1-2*varphi_w*frac).^0.5);
    scnd_term = 0;%g_e./(ne_0.*ce0).*sqrt(T_w./Te).*par;
    thrd_term = Mi0^2*((1-2.*varphi_w./Mi0^2).^0.5 - 1);
    half_dvarphi_dxi_squared = exp(varphi_w)-1 + scnd_term + thrd_term;
%     if half_dvarphi_dxi_squared<0
%        display(half_dvarphi_dxi_squared)
%        display(scnd_term)
%        display(thrd_term)
%     end
    dvarphi_dxi = (2*half_dvarphi_dxi_squared).^0.5;
    E_w = sqrt(ne_0.*Te/eps_0).*dvarphi_dxi;
end