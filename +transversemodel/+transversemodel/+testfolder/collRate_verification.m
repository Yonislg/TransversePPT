% test collision rate
import transversemodel.subfunctions.*;
Z=1;
ne_0 = 2*10^22;           % Electron density  [m^-3]
m_e = 9.1093837015e-31;
nb=ne_0;
n_i = ne_0/Z;
e = 1.602176634e-19;
Te = 2.4*e;               % Electron Temperature in joules (not eV!)
%vth_e = sqrt(8*Te/m_e/pi)
%u_xe = -5000;

nu_ei = collRate_ei(n_i,Z,Te) % search Markusic for typical values

eta = m_e*nu_ei/(e^2*ne_0)

%bulk_drop_N = u_xe*m_e*nb*(collRate_ei(nb,1,vth_e)) % N