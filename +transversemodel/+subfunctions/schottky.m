function ge = schottky(T_w, W, E_w, A_G)
    k = physconst('boltzmann');     % Bolzmann Constant [J/K]
    e = 1.602176634e-19;
    
    W = W*e; % ones(size(E_w))*
    delta_W = schottky_decrease(E_w);
    %delta_W = [0, 1.07, 1.57, 1.81, 2.01, 2.18]*e
    exponent = -(W - delta_W) /(k *T_w);
    ge = A_G *T_w.^2 / e .* exp(exponent); % including the term e depends on parameter A_G. Devide by e because to get flux
end


function delta_W = schottky_decrease(E_w)
    eps_0 = 8.854187817620389e-12; 
    e = 1.602176634e-19;
    pi_ac=3.14159265;
    delta_W = (e^3 * E_w / (4 * pi_ac * eps_0)).^0.5;
end