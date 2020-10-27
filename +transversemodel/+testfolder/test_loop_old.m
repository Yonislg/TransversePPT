%% Test loop 
% This script iterates over the relevant input parameters

close all
clear all
import transversemodel.subfunctions.*;
import  transversemodel.main.*;

%Physical constants
global k e eps_0 hbar m_e mu_0 u_ze imeq F_SEE F_FEE F_TEE;

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
T_iter=9;
n_iter = 20;

r_T_wka = [300 400 500 600 700];          % Anode temperture [k]
r_T_wkc = [500 600 700 800];          % Cathode temperture [k]
r_ne_0 = linspace(10^20,2*10^23, n_iter);           % Electron density  [m^-3]  2*10^22 3*10^22 4*10^22 
r_Te = linspace(e,3*e,T_iter);               % Electron Temperature in joules (not eV!)
r_phi_A =  [10 100 1000];          % Anode potential 
r_h = [0.01 0.02 0.05];               % Distance between electrodes
r_C_guess = [0 5 10 30];

iter = 200
%% Setting up tables
InputVarNames = {'Electron density', 'Electrode temperature [eV]','Anode potential' , 'Cathode temperature', 'Anode temperature', 'Electrode distance'};
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

x0=10;
y0=300;
width=670;
height=490;
%set(gcf,'position',[x0,y0,width,height])

%potential = figure;
%elec_vel = figure;
% figure(1)
% clf
% 
% title('Transverse Plasma Potential (Te=1 eV)')
% xlabel('Transverse position')
% ylabel('Potential [V]')
% hold on
% 
% figure(2)
% clf
% 
% title('Transverse Plasma Potential (Te=2 eV)')
% xlabel('Transverse position')
% ylabel('Potential [V]')
% hold on
% 
% figure(3)
% clf
% 
% title('Transverse Plasma Potential (Te=3 eV)')
% xlabel('Transverse position')
% ylabel('Potential [V]')
% hold on


figure(1)
clf
%  x0=620;
%  y0=300;
% set(gcf,'position',[x0,y0,width,height])
title('Electron Bulk Velocity')
xlabel('Electrode temperature [eV]')
ylabel('U_xe [m/s]')
 hold on
% 
figure(2)
clf
%  x0=620;
%  y0=300;
% set(gcf,'position',[x0,y0,width,height])
title('Electron Current')
xlabel('Electrode temperature [eV]')
ylabel('I [A/m^2]')
 hold on
% 

figure(3)
clf
title('Cathode Potentail Drop')
xlabel('Electrode temperature [eV]')
ylabel('Potential Difference [V]')
 hold on

figure(4)
clf
title('Bulk Potentail Drop')
xlabel('Electrode temperature [eV]')
ylabel('Potential Difference [V]')
 hold on
% figure(1)
% x0=10;
% y0=20;
% set(gcf,'position',[x0,y0,width,height])
% title('Electron Emissions')
% xlabel('Electrode temperature [eV]')
% ylabel('Total emissions [m^-2/s]')
% hold on
% 
% figure(3)
% x0=620;
% y0=20;
% set(gcf,'position',[x0,y0,width,height])
% %f4 = figure;
% title('Sheath Electric Fields')
% xlabel('Electrode temperature [eV]')
% ylabel('Electric field [v/m]')
% hold on



for p=1
    for j=1:n_iter
        if j==1
            style3='--';
        elseif j==2
            style3=':';
        end
        for i=1:T_iter
            for m= 2
                for n=3   
                    for v=3
    Te=r_Te(i);
    
    T_wka = r_T_wka(n); 
    T_wkc = r_T_wkc(m);
    ne_0 = r_ne_0(j);
    phi_A = r_phi_A(v);
    h = r_h(3);
    C_guess = r_C_guess(p);
    
    ui0 = sqrt(Te/m_i);      % Ion sheath boundary velocity
    

    plasma_properties = {Te, ne_0,n_n,Z, m_i};
    design_parameters = {T_wka, T_wkc, E_i, A_G, h, L, W, E_F};
    
    [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D,x,fx,exitflag,initial_state] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C,C_guess);

    if exitflag>0&&isreal(x)
        exitflag;
        ctr=ctr+1;
        InpuTable(ctr,:) = {ne_0, Te/e,phi_A, T_wkc, T_wka, h};
        IniTable(ctr,:) = array2table([initial_state phi_A-initial_state(2) phi_C-initial_state(1)]);
        OutpuTable(ctr,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D};
        ExiTable(ctr,:) = array2table([exitflag fx]);
        
        %plot potential drops
        
%         figure(i)
%         x = [0 0.25 2.75 3];
%         y = [phi_A OutpuTable.("Anode potential")(ctr) OutpuTable.("Cathode potential")(ctr) phi_C];
%         plot(x,y,style3,'DisplayName', sprintf('n = %.1e m^{-3}', ne_0),'color',cc(j,:))
%         xticks(x)
%         xticklabels({'A','B','D','C'})
        
        
    end
                    end
                end
            end
        end
        % plot all other variables
        % Electron bulk velocity
        figure(1)
        plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),-OutpuTable.('electron bulk velocity')(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', ne_0),'color',cc(j,:))
        
        figure(2)
        plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),-OutpuTable.('electron bulk velocity')(strt:ctr).*ne_0*e,style3,'DisplayName', sprintf('n = %.1e m^{-3}', ne_0),'color',cc(j,:))
        
        figure(3)
        plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.("Cathode potential")(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', ne_0),'color',cc(j,:))
        
        figure(4)
        plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.("Anode potential")(strt:ctr)-OutpuTable.("Cathode potential")(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', ne_0),'color',cc(j,:))
%         
%         %plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),IniTable.('electron bulk velocity')(strt:ctr),'--')
%         
        %Emission levels
%         if j==1
%             style1='b--';
%             style2='r:*';
%         elseif j==2
%             style1='g--';
%             style2='k:*';
%         end
%         figure(1)
%         %plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.('Cathode emissions')(strt:ctr),style1,InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.('Anode emissions')(strt:ctr),style2)  
%         plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.('Cathode emissions')(strt:ctr),style1)
%         plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.('Anode emissions')(strt:ctr),style2)  
%         plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),IniTable.('Cathode emissions')(strt:ctr))
%         plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),IniTable.('Anode emissions')(strt:ctr))  
        
%         % Electric fields
%         figure(3)
%         plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.('Cathode E-field')(strt:ctr),InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.('Anode E-field')(strt:ctr)) 
%         %plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),IniTable.('Cathode E-field')(strt:ctr),':',InpuTable.('Electrode temperature [eV]')(strt:ctr),IniTable.('Anode E-field')(strt:ctr),':') 

        strt = ctr+1;
    end
end

dens =cellfun(@(c) sprintf('%0.1e',c),num2cell(r_ne_0),'UniformOutput',false);

for i = 1:4
    figure(i)
    colormap(cc)
    hc= colorbar;
    set(hc,'Ticklabels',dens,'limit',[0 1],'Ticks',linspace(0,1,n_iter))
    set(hc.Label,'String','Density in m^{-3}')
    
end
%figure(1)
%legend(sprintf('n = %.0e m^{-3}, T = 1 ev', r_ne_0(1)) ,sprintf('n = %.0e m^{-3}, T = 2 ev', r_ne_0(1)), sprintf('n = %.0e m^{-3}, T = 3 ev', r_ne_0(1)),sprintf('n = %.0e m^{-3}, T = 1 ev', r_ne_0(2)),sprintf('n = %.0e m^{-3}, T = 2 ev', r_ne_0(2)), sprintf('n = %.0e m^{-3}, T = 3 ev', r_ne_0(2)))

% figure(1)
% %legend( {sprintf('n = %.0e m^{-3}', r_ne_0(1)),sprintf('n = %.0e m^{-3} (Init)', r_ne_0(1)),sprintf('n = %.0e m^{-3}', r_ne_0(2)),sprintf('n = %.0e m^{-3} (Init)', r_ne_0(2))},'Location','northeast');
% legend( {sprintf('n = %.0e m^{-3}', r_ne_0(1)),sprintf('n = %.0e m^{-3}', r_ne_0(2))},'Location','northeast');


% figure(1)
% legend({sprintf('Cathode n = %.0e m^{-3}', r_ne_0(1)),sprintf('Anode n = %.0e m^{-3}', r_ne_0(1)),sprintf('Cathode n = %.0e m^{-3} (init)', r_ne_0(1)),sprintf('Anode n = %.0e m^{-3} (init)', r_ne_0(1)),sprintf('Cathode n = %.0e m^{-3}', r_ne_0(2)),sprintf('Anode n = %.0e m^{-3}', r_ne_0(2)),sprintf('Cathode n = %.0e m^{-3} (init)', r_ne_0(2)),sprintf('Anode n = %.0e m^{-3} (init)', r_ne_0(2))},'Location','northwest');

% figure(3)
% 
% %legend( {sprintf('Cathode n = %.0e m^{-3}', r_ne_0(1)),sprintf('Anode n = %.0e m^{-3}', r_ne_0(1)), sprintf('Cathode n = %.0e m^{-3} (init)', r_ne_0(1)),sprintf('Anode n = %.0e m^{-3} (init)', r_ne_0(1)),sprintf('Cathode n = %.0e m^{-3}', r_ne_0(2)),sprintf('Anode n = %.0e m^{-3}', r_ne_0(2)),sprintf('Cathode n = %.0e m^{-3} (init)', r_ne_0(2)),sprintf('Anode n = %.0e m^{-3} (init)', r_ne_0(2))},'Location','northwest');
% legend( {sprintf('Cathode n = %.0e m^{-3}', r_ne_0(1)),sprintf('Anode n = %.0e m^{-3}', r_ne_0(1)),sprintf('Cathode n = %.0e m^{-3}', r_ne_0(2)),sprintf('Anode n = %.0e m^{-3}', r_ne_0(2))},'Location','northwest');

InpuTable(ctr+1:end,:)=[]
IniTable(ctr+1:end,:)=[]
OutpuTable(ctr+1:end,:)=[]
ExiTable(ctr+1:end,:)=[];

%% Make Matrices

OuTe = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),InpuTable.('Electrode temperature [eV]'));
CatDrop = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Cathode sheath potential drop'));
AnoDrop = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Anode sheath potential drop'));
CatEm = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Cathode emissions'));
AnoEm = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Anode emissions'));
CatEfield = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Cathode E-field'));
AnoEfield = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Anode E-field'));
Bulkvel = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('electron bulk velocity'));
CatPot = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Cathode potential'));
AnoPot = Tab2Mat(r_ne_0,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]'),OutpuTable.('Anode potential'));
ResMat = AnoPot-CatPot;



%% Colorising figures
dens =cellfun(@(c) sprintf('%0.1e',c),num2cell(r_ne_0),'UniformOutput',false);

for i = 1:4
    figure(i)
    colormap(cc)
    hc= colorbar;
    set(hc,'Ticklabels',dens,'limit',[0 1],'Ticks',linspace(0,1,n_iter))
    set(hc.Label,'String','Density in m^{-3}')
    
end


%% figures on resistivity
%Count the ammount of previous figures
h =  findobj('type','figure');
S = length(h);

% parameters
nb = InpuTable.('Electron density');
Te = InpuTable.('Electrode temperature [eV]');
eta_i = m_e*(collRate_ei(nb,1,Te*e))./nb/e^2;
eta_n =m_e*10^(-20)*sqrt(Te*e*m_e)/e^2;

% %figures
% figur(S+1)
% hold on
% for i = 1:length(T_iter)%n_iter
%     l=nb==r_Te(i)
%     plot(nb(l),eta_i(1:3:30),':',nb(1:3:30),eta_i(2:3:30),':',nb(1:3:30),eta_i(3:3:30),':')
%     plot(nb(1:3:30),eta_n(1:3:30),nb(1:3:30),eta_n(2:3:30),nb(1:3:30),eta_n(3:3:30))
% end


%legend('Te=1',"Te=2","Te=3",'Te=1',"Te=2","Te=3")
% figure(S+2)
%  R= m_e*(collRate_ei(nb,1,Te*e)).*nb.*-OutpuTable.('electron bulk velocity');
%  R2 = m_e*(nb.*sqrt(Te*e/m_e)*10^(-20)).*nb.*-OutpuTable.('electron bulk velocity');
% plot(nb(1:3:30),R(1:3:30),':',nb(1:3:30),R(2:3:30),':',nb(1:3:30),R(3:3:30),':')
% hold on
% plot(nb(1:3:30),R2(1:3:30),nb(1:3:30),R2(2:3:30),nb(1:3:30),R2(3:3:30))
% legend('Te=1',"Te=2","Te=3",'Te=1',"Te=2","Te=3")
% xlabel("Density")
% ylabel("Resistance")
% 
% figure(S+3)
% semilogy(nb(1:3:30),R(1:3:30),':',nb(1:3:30),R(2:3:30),':',nb(1:3:30),R(3:3:30),':',nb(1:3:30),R2(1:3:30),nb(1:3:30),R2(2:3:30),nb(1:3:30),R2(3:3:30))
% legend('Te=1',"Te=2","Te=3",'Te=1',"Te=2","Te=3")
% xlabel("Density")
% ylabel("Resistance")

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
