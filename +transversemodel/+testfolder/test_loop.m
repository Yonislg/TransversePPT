%% Test loop 
% This script iterates over the relevant input parameters

close all
clear all

import transversemodel.subfunctions.*;
import  transversemodel.main.*;
%Physical constants
global k e eps_0 hbar m_e mu_0 imeq F_SEE F_FEE F_TEE;


% To be added :  u_ze 

k = physconst('boltzmann');     % Bolzmann Constant [J/K]
e = 1.602176634e-19;
eps_0 = 8.854187817620389e-12; 
m_e = 9.1093837015e-31;         
m_i = 1.67262192369e-27;        % Ion mass, assumed to be equal to one proton
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
A_G = 80*10^4;          % Material constant for Schottkey equation

% parameter ranges:
T_iter= 25;
n_iter = 25;

r_T_wka = [300 400 500 600 700];                    % Anode temperture [k]
r_T_wkc = [500 600 700 800];                        % Cathode temperture [k]

LOGDENS = 0; % Set to 1 to plot density logarithmicaly. Set to 0 for linear

if LOGDENS
    r_nb = logspace(10^20,2*10^23, n_iter);           % Electron bulk density  [m^-3]  2*10^22 3*10^22 4*10^22 
else
    r_nb = linspace(10^20,2*10^23, n_iter);
end
r_Te = linspace(e,3*e,T_iter);                      % Electron Temperature in joules (not eV!)
r_phi_A =  [10 100 1000];                           % Anode potential 
r_h = [0.01 0.15 0.02 0.3 0.4 0.05];                             % Distance between electrodes
r_C_guess = [0 5 10 30];
u_ze = 10^4;             % Downstream (axial) flow velocity in m/s
By = 0.1;               % Magnetic field in Tesla


% Ionisation parameters 
a_iz = 0.5;     %ionisation degree
Z = 1;          %Ion charge number

iter = 2500

%% Setting up tables
InputVarNames = {'Electron density [m^{-3}]', 'Electrode temperature [eV]','Anode potential' , 'Cathode temperature', 'Anode temperature', 'Electrode distance'};
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
% 
% x0=10;
% y0=300;
% width=670;
% height=490;
% set(gcf,'position',[x0,y0,width,height])


for j = 1:n_iter
    for te = 1:T_iter
        Te=r_Te(te);
        T_wka = r_T_wka(3);
        T_wkc = r_T_wkc(2);
        nb = r_nb(j);       % electron bulk density
        %n_n = (1-a_iz)/a_iz*ne_0;
        phi_A = r_phi_A(3);
        h = r_h(6);
        C_guess = r_C_guess(1);
        
        plasma_properties = {Te, nb,a_iz,Z, m_i};
        design_parameters = {T_wka, T_wkc, E_i, A_G, h, L, W, E_F};
        
        [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx,exitflag,initial_state] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C, u_ze, By, C_guess);
        
        if exitflag>0&&isreal(x)
            exitflag;
            ctr=ctr+1;
            InpuTable(ctr,:) = {nb, Te/e,phi_A, T_wkc, T_wka, h};
            IniTable(ctr,:) = array2table([initial_state phi_A-initial_state(2) phi_C-initial_state(1)]);
            OutpuTable(ctr,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D};
            ExiTable(ctr,:) = array2table([exitflag fx]);
        end
    end
    strt = ctr+1;
end

InpuTable(ctr+1:end,:)=[];
IniTable(ctr+1:end,:)=[];
OutpuTable(ctr+1:end,:)=[];
ExiTable(ctr+1:end,:)=[];

%% Make Matrices

%matx = [r_nb,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]')]; 
label_1 = 'Electron density [m^{-3}]';
label_2 = 'Electrode temperature [eV]';

inmat1 = InpuTable.(label_1); 
inmat2 = InpuTable.(label_2);
invec1 = unique(inmat1);
invec2 = unique(inmat2);

NT = OrganizeFinds(inmat1,inmat2, InpuTable,OutpuTable);

%% Plot with invec2 on the x-axis and color coding representing invec1

% create color scheme
cc=jet(n_iter);



figure(1)
%plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),-OutpuTable.('electron bulk velocity')(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))
    hbv = plot(invec2,-NT.Bulkvel');
    set(hbv, {'color'},num2cell(cc,2))
    title('Electron Bulk Velocity')
    xlabel(label_2)
    ylabel('U_xe [m/s]')

figure(2)
%plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),-OutpuTable.('electron bulk velocity')(strt:ctr).*nb*e,style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))
    hcu = plot(invec2,(-NT.Bulkvel.*r_nb'*e)');
    set(hcu, {'color'},num2cell(cc,2))
    title('Electron Current')
    xlabel(label_2)
    ylabel('I [A/m^2]')

figure(3)
%plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.("Cathode potential")(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))
    hcat = plot(invec2,NT.CatPot');
    set(hcat, {'color'},num2cell(cc,2))
    title('Cathode Potential Drop')
    xlabel(label_2)
    ylabel('Potential Difference [V]')

figure(4)
    hres = plot(invec2,NT.ResMat');
    set(hres, {'color'},num2cell(cc,2))
    title('Bulk Potentail Drop')
    xlabel(label_2)
    ylabel('Potential Difference [V]')
    
%plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.("Anode potential")(strt:ctr)-OutpuTable.("Cathode potential")(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))



% set ticks labaels for for invec1
if LOGDENS % in the case invec1= density and logarithmically plotted
    tiLabls1 = cellfun(@(c) sprintf('%0.1e',c),num2cell(logspace(invec1(1),invec1(end),11)),'UniformOutput',false)
else
    tiLabls1 = cellfun(@(c) sprintf('%0.1e',c),num2cell(linspace(invec1(1),invec1(end),11)),'UniformOutput',false);
end


for i = 1:4
    figure(i)
    colormap(cc)
    hc= colorbar;
    set(hc,'Ticklabels',tiLabls1,'limit',[0 1],'Ticks',linspace(0,1,11))
    set(hc.Label,'String','Density in m^{-3}')
end


%% Plot with invec1 on the x-axis and color coding representing invec2


h =  findobj('type','figure');
S = 4%length(h);
c2=jet(length(unique(inmat2)));
    figure(S+1)
    hrm = plot(invec1,NT.ResMat);
    set(hrm, {'color'},num2cell(c2,2))
    figure(S+2)
    hcp = plot(invec1,NT.CatPot);
    set(hcp, {'color'},num2cell(c2,2))
    figure(S+3)
    hce = plot(invec1,-NT.Bulkvel);
    set(hce, {'color'},num2cell(c2,2))
    figure(S+4)
    hce = plot(invec1,-NT.Bulkvel.*r_nb'*e);
    set(hce, {'color'},num2cell(c2,2))


figure(S+1)
title('Bulk resistance')
xlabel(label_1)
ylabel('Potentail Drop [V]')

figure(S+2)
title('Cahtode Potential Drop')
xlabel(label_1)
ylabel('Potentail Drop [V]')
% 
figure(S+3)
%title('Cathode Emissions')
title('Electron Bulk Velocity')
xlabel(label_1)
ylabel('U_xe [m/s]')
%ylabel('Electron flux [m^{-3}/s]')

figure(S+4)
title('Electron Current')
xlabel(label_1)
ylabel('I [A/m^2]')



% Coroize temperaturs

%Temps = cellfun(@(c) sprintf('%0.1e',c),num2cell([r_Te(1)/e r_Te(T_iter/10:T_iter/10:T_iter)/e]),'UniformOutput',false); % make robost to change

tiLabls2 = cellfun(@(c) sprintf('%0.1e',c),num2cell(linspace(invec2(1),invec2(end),11)),'UniformOutput',false); % 

for i = S+1:S+4
    figure(i)
    colormap(c2)
    hc = colorbar;
    set(hc,'Ticklabels',tiLabls2,'limit',[0 1],'Ticks',linspace(0,1,11))
    set(hc.Label,'String','Temperature in eV')
    
end


%% 3D plot to check
S = 8;
[X,Y] = meshgrid(invec1, invec2);
C = X.*(Y*e/m_i).^0.5;

    figure(S+1)
    surf(X,Y,NT.ResMat',C);

    figure(S+2)
    surf(X,Y,NT.CatPot',C);

    figure(S+3)
    surf(X,Y,-NT.Bulkvel',C);

    figure(S+4)
    surf(X,Y,(-NT.Bulkvel.*r_nb'*e)',C);




figure(S+1)
title('Bulk resistance')
xlabel('Density [m^{-3}]')
ylabel('Electron temperature [eV]')

figure(S+2)
title('Cahtode Potential Drop')
xlabel('Density [m^{-3}]')
ylabel('Electron temperature [eV]')
% 
figure(S+3)
title('Electron Bulk Velocity')
xlabel('Density [m^{-3}]')
ylabel('Electron temperature [eV]')

figure(S+4)
title('Electron Current')
xlabel('Density [m^{-3}]')
ylabel('Electron temperature [eV]')

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
%% Reorganising outputs into a structure

function OutStruct = OrganizeFinds(inmat1,inmat2, InpuTable,OutpuTable)
    OutStruct.OuTe = Tab2Mat(inmat1, inmat2, InpuTable.('Electrode temperature [eV]'));
    OutStruct.CatDrop = Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode sheath potential drop'));
    OutStruct.AnoDrop = Tab2Mat(inmat1, inmat2, OutpuTable.('Anode sheath potential drop'));
    OutStruct.CatEm = Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode emissions'));
    OutStruct.AnoEm = Tab2Mat(inmat1, inmat2, OutpuTable.('Anode emissions'));
    OutStruct.CatEfield = Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode E-field'));
    OutStruct.AnoEfield = Tab2Mat(inmat1, inmat2, OutpuTable.('Anode E-field'));
    OutStruct.Bulkvel =  Tab2Mat(inmat1, inmat2, OutpuTable.('electron bulk velocity'));
    OutStruct.CatPot  =  Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode potential'));
    OutStruct.AnoPot  =  Tab2Mat(inmat1, inmat2, OutpuTable.('Anode potential'));
    OutStruct.ResMat  =  OutStruct.AnoPot-OutStruct.CatPot;
end

%% Matrix Maker
function OutMat = Tab2Mat(InCol1,InCol2,ProCol)
    vec1=unique(InCol1);
    vec2=unique(InCol2);
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
