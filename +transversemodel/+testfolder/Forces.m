%% This script examines the balance of forces
% At this point By is taken as an unkown, later tests will take By as a
% function of uxe
import transversemodel.subfunctions.*;

 % equation to be evaluated:
%jxe_v2(7) = -e*nb*(phi_A-x(2)*Te/e-(phi_C-x(1)*Te/e))/h +e*nb*u_ze*By - x(7)*m_e*nb*(collRate_ei(ni_b,Z,Te)+n_n*vth_e*10^(-20));

%% Start 
% Physiscs constsants
e = 1.602176634e-19;            % Electron charge
% parameters
Te = 3*e;
Z = 2;
a_iz = 1;
V= logspace(1,3,5);%1300;
u_ze = logspace(2,4,5);%10^4;
By = linspace(0.01,0.1,5);%logspace(-2,-1,20);
h=0.05;
nb=logspace(22,24,5);
L=0.001;


uxe = force_balance(V,nb,Z, a_iz,h,Te,By,u_ze)

plot(By,-uxe.*nb*e)
%I = -nb*e*force_balance2(V,nb,Z, a_iz,h,Te,L,u_ze)

%% rewriten as a function for uxe:
function uxe = force_balance(V,nb,Z, a_iz,h,Te,By,u_ze)
    import transversemodel.subfunctions.*;
    e = 1.602176634e-19;   
    ni_b= nb/Z;
    n_n = (1-a_iz)/a_iz*nb;
    m_e = 9.1093837015e-31;         % Electron mass
    vth_e = sqrt(Te/m_e);
    
    uxe = (-e*V/h +e*u_ze.*By)./(m_e*(collRate_ei(ni_b,Z,Te)+n_n*vth_e*10^(-20)));
end

function uxe = force_balance2(V,nb,Z, a_iz,h,Te,L,u_ze)
    import transversemodel.subfunctions.*;     
    mu_0 = 1.2566370614359173e-06;
    e = 1.602176634e-19;   
    ni_b= nb/Z;
    n_n = (1-a_iz)/a_iz*nb;
    m_e = 9.1093837015e-31;         % Electron mass
    vth_e = sqrt(Te/m_e);
    
    uxe = e*V/h/(e^2*u_ze*nb*mu_0*L -m_e*(collRate_ei(ni_b,Z,Te)+n_n*vth_e*10^(-20)));
end