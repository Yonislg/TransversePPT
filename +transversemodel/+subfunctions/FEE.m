%% Fowler-Nordheim
% Remove W from globals at later stage
function ge_FEE = FEE(W,E_F, E_W,CorF)
if nargin==3
    CorF=1;
end

e = 1.602176634e-19;
hbar = 1.0546e-34;
m_e = 9.1093837015e-31;
%W=W*e;
%E_F=E_F*e;

    a = e^2/((4*pi)^2*hbar);     % 1.5414e-06 A eV V^-2
    a=a/e;                       % Conversion to m^-3s^-1 eV V^-2    
    b = 4*sqrt(2*m_e)*e^(3/2)/(3*e*hbar);  %  6.8307 eV^-3/2 V nm^-1
    PF = 4*W/(W+E_F)*sqrt(E_F/W); %from Fridman and R. Forbes (Sommerfeld model)
    corF1 = PF;% 
    corF2 = v(E_W,W); % Correction factor, dependent on Delta_W/W , now asumed to be 1% CorF;%corF;
    ge_FEE = corF1.*E_W.^2*a/W.*exp(-b*W^(3/2)*corF2./E_W);

end

function corF = v(E_W,W)
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
    c =sqrt(e/(4*pi*eps_0));
    y = c*sqrt(E_W)/W;
    if y==0
        corF = 1;
    else
        corF = 1- y.^2+(1/3).*y.^2.*log(y);
    end
end
