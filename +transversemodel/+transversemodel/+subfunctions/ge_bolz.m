% Calculate sheath electron current using bolzmann
% flux = 1/4*v_e*exp(-eV/Te)

function I = ge_bolz(ne, Te, V)
    global e m_e;
    %ne * (Te / (2 * pi * m_e)).^ 0.5
    
    I = ne * (Te / (2 * pi * m_e)) ^ 0.5 * exp(-e * V / Te);
end