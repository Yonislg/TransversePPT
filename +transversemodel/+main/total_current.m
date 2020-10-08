% This equation calculates the current in a 0D system using 7 eq's provided
% by mario

function jxe_v2 = total_current(x, plasma_properties,design_parameters, phi_A, phi_C)
    import transversemodel.subfunctions.*;
    global m_i mu_0 m_e e u_ze imeq F_SEE F_FEE;


    [Te, ne_0, ui0] = deal(plasma_properties{:});
    nb=ne_0;
    [T_wka, T_wkc, E_i, A_G, h, L, W, E_Fin] = deal(design_parameters{:});
if F_SEE
    ge_SEE = SEE(ne_0*sqrt(Te/m_i), E_i, W);
else 
    ge_SEE = 0;
end

if F_FEE
    E_F= E_Fin;
else 
    E_F = 0;
end
    vth_e = sqrt(Te/m_e);% 1D RMS Thermal velocity of electrons
    % Cathode electric field
    jxe_v2(1) = wall_e_field(T_wkc, x(1),x(3),ne_0, Te) - x(5); 
    %Cathode lectron emissions
    jxe_v2(2) = - x(3)  + schottky(T_wkc, W, x(5), A_G)+ ge_SEE + FEE(W,E_F, x(5));
    % Average of fluxes
    if imeq
        jxe_v2(3) = ((x(4) - ge_bolz(ne_0, Te, -x(2)*Te/e)  - x(3) + ge_bolz(ne_0, Te, -x(1)*Te/e))/2)/ne_0 - x(7);
    else
        jxe_v2(3) = (x(4) - ge_bolz(ne_0, Te, -x(2)*Te/e)  - x(3) + ge_bolz(ne_0, Te, -x(1)*Te/e))/2/ne_0;
    end
    % Continuity equation
    jxe_v2(4) = -2*sqrt(Te/m_i) + (ge_bolz(ne_0, Te, -x(2)*Te/e)+ ge_bolz(ne_0, Te, -x(1)*Te/e)- x(4)- x(3))/ne_0;
    % Anode electron emissions
    jxe_v2(5) =  - x(4)+ schottky(T_wka, W, x(6), A_G) + ge_SEE + FEE(W,E_F, x(6)); 
    % Anode electric field
    jxe_v2(6) = wall_e_field(T_wka, x(2),x(4),ne_0, Te) -x(6);  
    % Momentum equation
    if imeq
        jxe_v2(7) = -e*nb*(phi_A-x(2)*Te/e-(phi_C-x(1)*Te/e))/h - x(7)*m_e*nb*(collRate_ei(nb,1,Te)+nb*vth_e*10^(-20));%+ u_ze*(e*ne_0*L*x(7)*mu_0); % Bulk current
    end
end