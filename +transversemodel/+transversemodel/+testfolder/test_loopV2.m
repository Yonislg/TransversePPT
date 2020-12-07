%% Test loop 
% This script iterates over the relevant input parameters

close all
clear all

import transversemodel.subfunctions.*;
import  transversemodel.main.*;
%Physical constants
global k e eps_0 hbar m_e m_i mu_0 u_ze imeq F_SEE F_FEE F_TEE;

k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;
m_i = 1.67262192369e-27;
mu_0 = 1.2566370614359173e-06;
hbar = 1.0546e-34;

% Flags
imeq=1; %Include momentum balance equation 
if imeq
    fprintf('Momentum equation is included');
else
    fprintf('Momentum equation is not included');
end

% Flags
F_SEE=1; %Include SEE 
if F_SEE
    fprintf('\nSecondary electron emmissions are included');
else
    fprintf('\nSecondary electron emmissions are not included');
end

F_TEE=0; %Include SEE 
if F_TEE
    fprintf('\nThermal electron emmissions are included');
else
    fprintf('\nThermal electron emmissions are not included');
end

F_FEE=0; %Include FEE 
if F_FEE
    fprintf('\nField electron emmissions are included');
else
    fprintf('\nField electron emmissions are not included');
end

%fixed 
W= 4.4; % Workfunction of copper (4.4) in eV
E_i = 13.6;             % Ionization energy of hydrogen in electronvolt [eV]
E_F = 7;                % Fermi energy [eV]
phi_C = 0;           % Cathode potenital
L = 0.08;      
A_G = 80*10^4;

% parameter ranges:
T_iter=3;
n_iter = 20;
iziter = 20;

r_T_wka = [300 400 500 600 700];          % Anode temperture [k]
r_T_wkc = [500 600 700 800];          % Cathode temperture [k]
r_ne_0 = logspace(7,22, n_iter);           % Electron density  [m^-3]  2*10^22 3*10^22 4*10^22 
r_Te = linspace(e,3*e,T_iter);               % Electron Temperature in joules (not eV!)
r_phi_A =  [10 100 1000];          % Anode potential 
r_h = [0.01 0.02 0.05];               % Distance between electrodes
r_C_guess = [0 5 10 30];
% Ionisation parameters 
r_a_iz = linspace(0.001,0.2,iziter)
%a_iz = 0.5;     %ionisation degree
Z = 1;          %Ion charge number

iter = 2500
%% Setting up tables
InputVarNames = {'Electron density', 'Electrode temperature [eV]','Ionization degree \alpha_{iz}','Anode potential' , 'Cathode temperature', 'Anode temperature', 'Electrode distance'};
varTypes=repmat({'double'},1,length(InputVarNames));
sz= [iter length(InputVarNames)];
InpuTable = table('Size', sz, 'VariableTypes',varTypes, 'VariableNames', InputVarNames);

OutputVarNames = {'Cathode sheath potential drop','Anode sheath potential drop', 'Cathode emissions','Anode emissions', 'Cathode E-field','Anode E-field','electron bulk velocity','Anode potential', 'Cathode potential'};
varTypes2=repmat({'double'},1,length(OutputVarNames));
sz2= [iter length(OutputVarNames)];
IniTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);
OutpuTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);

ExitVarNames = {'ExitFlag','Cathode sheath potential drop','Anode sheath potential drop', 'Cathode emissions','Anode emissions', 'Cathode E-field','Anode E-field','electron bulk velocity'};
varTypes3=['int8' repmat({'double'},1,length(ExitVarNames)-1)];
sz3= [iter length(ExitVarNames)];

ExiTable=table('Size', sz3, 'VariableTypes',varTypes3, 'VariableNames', ExitVarNames);
%%
ctr = 0;
strt = 1;

% create color scheme
cc=jet(n_iter);

% 
% figure(1)
% clf
% %  x0=620;
% %  y0=300;
% % set(gcf,'position',[x0,y0,width,height])
% title('Electron Bulk Velocity')
% xlabel('Electron temperature [eV]')
% ylabel('U_xe [m/s]')
%  hold on
% % 
% figure(2)
% clf
% %  x0=620;
% %  y0=300;
% % set(gcf,'position',[x0,y0,width,height])
% title('Electron Current')
% xlabel('Electron temperature [eV]')
% ylabel('I [A/m^2]')
%  hold on
% % 
% 
% figure(3)
% clf
% title('Cathode Potentail Drop')
% xlabel('Electron temperature [eV]')
% ylabel('Potential Difference [V]')
%  hold on
% 
% figure(4)
% clf
% title('Bulk Potentail Drop')
% xlabel('Electron temperature [eV]')
% ylabel('Potential Difference [V]')
% hold on


for p=1
    for j=1:n_iter
        if j==1
            style3='--';
        elseif j==2
            style3=':';
        end
        for i=3
            for m= 2
                for n=3   
                    for v=3
                        for iz = 1:iziter
    Te=r_Te(i);
    
    T_wka = r_T_wka(n); 
    T_wkc = r_T_wkc(m);
    ne_0 = r_ne_0(j);
    a_iz = r_a_iz(iz);
    n_n = (1-a_iz)/a_iz*ne_0;
    phi_A = r_phi_A(v);
    h = r_h(3);
    C_guess = r_C_guess(p);
    
    %ui0 = sqrt(Te/m_i);      % Ion sheath boundary velocity
    

    plasma_properties = {Te, ne_0,n_n,Z};
    design_parameters = {T_wka, T_wkc, E_i, A_G, h, L, W, E_F};
    
    [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx,exitflag,initial_state] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C,C_guess);

    if exitflag>0&&isreal(x)
        exitflag;
        ctr=ctr+1;
        InpuTable(ctr,:) = {ne_0, Te/e, a_iz, phi_A, T_wkc, T_wka, h};
        IniTable(ctr,:) = array2table([initial_state phi_A-initial_state(2) phi_C-initial_state(1)]);
        OutpuTable(ctr,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D};
        ExiTable(ctr,:) = array2table([exitflag fx]);
               
        
    end
                        end
                    end
                end
            end
        end
        strt = ctr+1;
    end
end


InpuTable(ctr+1:end,:)=[]
IniTable(ctr+1:end,:)=[]
OutpuTable(ctr+1:end,:)=[]
ExiTable(ctr+1:end,:)=[];

%% Make Matrices

OuTe = Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),InpuTable.('Electrode temperature [eV]'));
CatDrop = Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Cathode sheath potential drop'));
AnoDrop = Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Anode sheath potential drop'));
CatEm = Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Cathode emissions'));
AnoEm = Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Anode emissions'));
CatEfield = Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Cathode E-field'));
AnoEfield = Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Anode E-field'));
Bulkvel =  Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('electron bulk velocity'));
CatPot  =  Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Cathode potential'));
AnoPot  =  Tab2Mat(r_ne_0,r_a_iz,InpuTable.('Electron density'),InpuTable.('Ionization degree \alpha_{iz}'),OutpuTable.('Anode potential'));
ResMat  =  AnoPot-CatPot;

% h =  findobj('type','figure');
% S = 4;%= length(h);
S=0;
c2=jet(iziter);
    figure(S+1)
    hrm = plot(r_ne_0,ResMat);
    set(hrm, {'color'},num2cell(c2,2))
    figure(S+2)
    hcp = plot(r_ne_0,CatPot);
    set(hcp, {'color'},num2cell(c2,2))
%     figure(S+3)
%     hce = plot(r_ne_0,CatEm);
%     set(hce, {'color'},num2cell(c2,2))

%
figure(S+1)
title('Bulk resistance')
xlabel('Density [m^{-3}]')
ylabel('Potentail Drop [V]')

figure(S+2)
title('Cahtode Potential Drop')
xlabel('Density [m^{-3}]')
ylabel('Potentail Drop [V]')

% figure(S+3)
% title('Cathode Emissions')
% xlabel('Density [m^{-3}]')
% ylabel('Electron flux [m^{-3}/s]')



%% Plot for Temps
Temps = cellfun(@(c) sprintf('%0.1e',c),num2cell(r_a_iz(iziter/10:iziter/10:iziter)),'UniformOutput',false); % make robost to change

for i = 1:2
    figure(i)
    colormap(c2)
    hc = colorbar;
    set(hc,'Ticklabels',Temps,'limit',[0 1],'Ticks',linspace(0,1,10))
    set(hc.Label,'String','\alpha_iz')%'Temperature in eV')
    
end

% %% Colorising figures
% dens =cellfun(@(c) sprintf('%0.1e',c),num2cell([r_ne_0(1) r_ne_0(n_iter/10:n_iter/10:n_iter)]),'UniformOutput',false);
% 
% for i = 1:4
%     figure(i)
%     colormap(cc)
%     hc = colorbar;
%     set(hc,'Ticklabels',dens,'limit',[0 1],'Ticks',linspace(0,1,11))
%     set(hc.Label,'String','Density in m^{-3}')
%     
% end
% 
% 
% %% figures on resistivity
% %Count the ammount of previous figures
% h =  findobj('type','figure');
% S = length(h);
% 
% % parameters
% nb = InpuTable.('Electron density');
% Te = InpuTable.('Electrode temperature [eV]');
% eta_i = m_e*(collRate_ei(nb,1,Te*e))./nb/e^2;
% eta_n =m_e*10^(-20)*sqrt(Te*e*m_e)/e^2;

% %figures
% figur(S+1)
% hold on
% for i = 1:length(T_iter)%n_iter
%     l=nb==r_Te(i)
%     plot(nb(l),eta_i(1:3:30),':',nb(1:3:30),eta_i(2:3:30),':',nb(1:3:30),eta_i(3:3:30),':')
%     plot(nb(1:3:30),eta_n(1:3:30),nb(1:3:30),eta_n(2:3:30),nb(1:3:30),eta_n(3:3:30))
% end


%% Matrix Maker
function OutMat = Tab2Mat(vec1,vec2,InCol1,InCol2,ProCol)
    Outmat1 = zeros(length(vec1),length(vec2));
    for i  = 1:length(vec1)
        for j = 1:length(vec2)
            if isempty(ProCol((InCol1==vec1(i))&(InCol2==vec2(j))))
                Outmat1(i,j) = NaN;
            else
                Outmat1(i,j) = ProCol((InCol1==vec1(i))&(InCol2==vec2(j)));
            end
        end
    end
    OutMat = Outmat1;
end
