function ge = schottky(T_w, W, E_w, A_G,fix)
    k = physconst('boltzmann');     % Bolzmann Constant [J/K]
    e = 1.602176634e-19;
    if exist('fix','var') 
        if fix>=0
        E_w = fix;9.3548565e+07;%2.74641711143e+07;%9.425091168*10^7;%9.4234e+07;%9.5548e+07;%9.42509*10^7;
        %fix = 0;
        end
    end
    
    %if fix
        
%     else
%         eps_0 = 8.854187817620389e-12; 
%         Te = 2*e;
%         ne_0 = 5*10^21;
%         E_w = sqrt(ne_0.*Te/eps_0)*E_w;
    %end
    W = W*e; % ones(size(E_w))*
    
    
    delta_W = schottky_decrease(E_w);
    
    
    %delta_W = [0, 1.07, 1.57, 1.81, 2.01, 2.18]*e
    exponent = -(W - delta_W)./(k .*T_w);
    ge = A_G .*T_w.^2 / e .* exp(exponent); % including the term e depends on parameter A_G. Devide by e because to get flux
end


function delta_W = schottky_decrease(E_w)
    eps_0 = 8.854187817620389e-12; 
    e = 1.602176634e-19;
    pi_ac=3.14159265;
    delta_W = (e^3 * E_w / (4 * pi_ac * eps_0)).^0.5;
end