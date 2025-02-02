function nu = collRate_ei(n_i,Z,T_e) 
    % n is the ion density!!!
    m_e = 9.1093837015e-31; % is correct
 
    e = 1.602176634e-19;
    eps_0 = 8.854187817620389e-12; 
    
    m=m_e;
    ne = n_i*Z;
    
    lnlambda= 23 - log(ne.^(1/2)*Z.*(T_e/e).^(-3/2))/log(10);%10;% Using 23 - ln(ne^1/2*Zi*Te^-3/2)

    %nu = n*Z^2*e^4*lnlambda/(4*pi*eps_0^2*m^2*v^3); Gold ston
    
    % From coulomb comparison, T_e in joules
    nu = sqrt(2).*n_i*Z^2*e^4.*lnlambda./(12*pi^(3/2)*eps_0^2*m^0.5.*T_e.^(3/2));
    
    %    1.8*10^5*(n*10^(-20))/T^(3/2);
end