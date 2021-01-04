%% Calulcates current density and sheath potential drop for transverse slabs of the plasma

function jxe_v2 = total_current(x, plasma_properties,design_parameters, phi_A, phi_C,u_ze, By)
    %% prepping script
    import transversemodel.subfunctions.*;
    global  imeq F_SEE F_FEE;
    % to be adde as input: u_ze, ion mass number.
    % Physiscs constsants
    e = 1.602176634e-19;            % Electron charge
    m_e = 9.1093837015e-31;         % Electron mass
    
    %loading from code
    [Te, ne_0, nb, n_n, Z, m_i] = deal(plasma_properties{:});
    % bulk density, probably double of sheath edge density
    [T_wka, T_wkc, E_i, A_G, h, L, W, E_Fin] = deal(design_parameters{:});
    ni_0= ne_0/Z;        % Ion density
    ni_b = nb/Z;
    % Check if SEE is to be included
if F_SEE
    ge_SEE = SEE(ni_0*sqrt(Te/m_i), E_i, W);
else 
    ge_SEE = 0;
end
    % Check if FEE is to be included

if F_FEE
    E_F= E_Fin;
else 
    E_F = 0;
end
    vth_e = sqrt(Te/m_e);% 1D RMS Thermal velocity of electrons
    
    %% the system of eqautions
    % Cathode electric field
    jxe_v2(1) = (wall_e_field(T_wkc, x(1),x(3),ne_0, Te) - x(5)); 
    %Cathode lectron emissions
    jxe_v2(2) = - x(3)  + schottky(T_wkc, W, x(5), A_G)+ ge_SEE + FEE(W,E_F, x(5));
    % Average of fluxes, when imeq is set two 0 the electron current is
    % ignored
    if imeq
        jxe_v2(3) = ((x(4) - ge_bolz(ne_0, Te, -x(2)*Te/e)  - x(3) + ge_bolz(ne_0, Te, -x(1)*Te/e))/2)/nb - x(7);
    else
        jxe_v2(3) = (x(4) - ge_bolz(ne_0, Te, -x(2)*Te/e)  - x(3) + ge_bolz(ne_0, Te, -x(1)*Te/e))/2/nb;
    end
    % Continuity equation
    jxe_v2(4) = -2*sqrt(Te/m_i) + (ge_bolz(ne_0, Te, -x(2)*Te/e)+ ge_bolz(ne_0, Te, -x(1)*Te/e)- x(4)- x(3))/ne_0;
    % Anode electron emissions
    jxe_v2(5) =  - x(4)+ schottky(T_wka, W, x(6), A_G) + FEE(W,E_F, x(6)) + ge_SEE ;
    % Anode electric field
    jxe_v2(6) = wall_e_field(T_wka, x(2),x(4),ne_0, Te) -x(6);  
    % Momentum equation
    if imeq
        jxe_v2(7) = -e*nb*(phi_A-x(2)*Te/e-(phi_C-x(1)*Te/e))/h +e*nb*u_ze*By - x(7)*m_e*nb*(collRate_ei(ni_b,Z,Te)+n_n*vth_e*10^(-20));%+ u_ze*(e*nb*L*x(7)*mu_0); % Bulk current
    end
end