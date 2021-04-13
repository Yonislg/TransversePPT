%% Test loop 
% This script iterates over the relevant input parameters

close all
clear all

import transversemodel.subfunctions.*;
import  transversemodel.main.*;
%Physical constants
global k e eps_0 hbar m_e mu_0 imeq;


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

% 
% % Flags
% F_SEE=1; %Include SEE 
% if F_SEE
%     fprintf('\nSecondary electron emmissions are included');
%     SEE_inc = 'yes';
% else
%     fprintf('\nSecondary electron emmissions are not included');
%     SEE_inc = 'a posteriori';
% end
% 
% F_TEE=0; %Include SEE 
% if F_TEE
%     fprintf('\nThermal electron emmissions are included');
%     TEE_inc = 'yes';
% else
%     fprintf('\nThermal electron emmissions are not included');
%     TEE_inc = 'a posteriori';
% end
% 
% F_FEE=0; %Include FEE 
% if F_FEE
%     fprintf('\nField electron emmissions are included');
%     FEE_inc = 'yes';
% else
%     fprintf('\nField electron emmissions are not included');
%     FEE_inc = 'a posteriori';
% end

%fixed 
W= 4.4; % Workfunction of copper (4.4) in eV
E_i = 13.6;             % Ionization energy of hydrogen in electronvolt [eV]
E_F = 7;                % Fermi energy [eV]
phi_C = 0;           % Cathode potenital
L = 0.08;      
A_G = 80*10^4;          % Material constant for Schottkey equation

% parameter ranges:
T_iter= 3;
n_iter = 3;

r_T_wka = [300 400 500 600 700];                    % Anode temperture [k]
r_T_wkc = linspace(3490.7122,3490.7123,11);
%[800 950 960 961 961.5 962];%1000 2000 3000 3400]%linspace(2000,3500,25); %[500 600 700 800 2500];                        % Cathode temperture [k]

LOGDENS = 0; % Set to 1 to plot density logarithmicaly. Set to 0 for linear

if LOGDENS
    r_nb = logspace(22,23, n_iter);           % Electron bulk density  [m^-3]  2*10^22 3*10^22 4*10^22 
else
    r_nb = linspace(10^22, 10^24, n_iter);
end
r_Te = linspace(e,3*e,T_iter);                      % Electron Temperature in joules (not eV!)
r_phi_A =  logspace(1,3,T_iter);%[10 100 1000];                           % Anode potential 
r_h = linspace(0.01, 0.05, n_iter);%[0.01 0.15 0.02 0.3 0.4 0.05];                             % Distance between electrodes
%r_B = linspace(0,0.8,n_iter);
r_C_guess = [0 1 2 5 10 30];
u_ze = 10^4;             % Downstream (axial) flow velocity in m/s

%By = 0.7;               % Magnetic field in Tesla

d_J = 0%0.01;

% Ionisation parameters 
a_iz = 1;     %ionisation degree
Z = 1;          %Ion charge number

iter = 60
% Iteration limit for outer loop on E_W
OutLim = 20;
epres = 10^-12;         % required precision for E_w

%A_G1 =A_G*F_TEE;%      % Material Constant for shottkey equation
%% Setting up tables
 InputVarNames = {'Electron density [m^{-3}]', 'Electrode temperature [eV]','Anode potential' , 'Cathode temperature', 'Anode temperature', 'Electrode distance','Slug Thickness'};
 varTypes=repmat({'double'},1,length(InputVarNames));
 sz= [iter length(InputVarNames)];
 InpuTable = table('Size', sz, 'VariableTypes',varTypes, 'VariableNames', InputVarNames);
 
OutputVarNames = {'Cathode sheath potential drop','Anode sheath potential drop', 'Cathode emissions','Anode emissions', 'Cathode E-field','Anode E-field','electron bulk velocity','Anode potential', 'Cathode potential','Magnetic Field'};
varTypes2=repmat({'double'},1,length(OutputVarNames));
sz2= [iter length(OutputVarNames)];
linit = length(OutputVarNames)-1;
sz4= [iter linit];
varTypes4=repmat({'double'},1,linit);
IniTable = table('Size', sz4, 'VariableTypes',varTypes4, 'VariableNames', OutputVarNames(1:linit));
 OutpuTable = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);

ExitVarNames = {'ExitFlag', 'Cathode sheath potential drop','Anode sheath potential drop', 'Cathode emissions','Anode emissions', 'Cathode E-field','Anode E-field','electron bulk velocity'};
varTypes3=['int8' repmat({'double'},1,length(ExitVarNames)-1)];
sz3= [iter length(ExitVarNames)];
ExiTable=table('Size', sz3, 'VariableTypes',varTypes3, 'VariableNames', ExitVarNames);

%% Setting up table
InputVarNames = {'n_e [m^{-3}]', 'T_e [eV]','Phi_A [V]' , 'T_{w,A} [K]', 'T_{w,c} [K]', 'h [m]','\delta [m]','SEE','TEE','FEE'};
varTypes=[repmat({'double'},1,length(InputVarNames)-3) repmat({'string'},1,3)];
sz= [iter length(InputVarNames)];
Input = table('Size', sz,'VariableTypes',varTypes, 'VariableNames', InputVarNames);

OutputVarNames = {'U_{CD} [V]','U_{AB} [V]', 'g*_C [m^2/s]','g*_A [m^2/s]', 'E_{w,c} [V/m]','E_{w,A} [V/m]','u_ze [m/s]','Phi_B', 'Phi_D','B_y [T]'};
varTypes2=repmat({'double'},1,length(OutputVarNames));
sz2= [iter length(OutputVarNames)];
Output = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', OutputVarNames);

PotVarNames = {'U_{CD} [V]','U_{AB} [V]','E_{w,c} [V/m]','E_{w,A} [V/m]','Phi_B [V]', 'Phi_D [V]','E_x [V/m]'};
varTypes2 = repmat({'double'},1,length(PotVarNames));
sz2= [iter length(PotVarNames)];
Potential = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', PotVarNames);

EmissioNames = {'g_C [a/m^2]','g*_C [a/m^2]', 'SEE_C [a/m^2]', 'TEE_C [a/m^2]','g_A [a/m^2]','g*_A [a/m^2]', 'SEE_A [a/m^2]', 'TEE_A [a/m^2]', 'J [a/m^2]'};
varTypes3=repmat({'double'},1,length(EmissioNames));
sz3= [iter length(EmissioNames)];
Emissions = table('Size', sz3, 'VariableTypes',varTypes3, 'VariableNames', EmissioNames);

ErrorNames = {'E_{w,c} [V/m]','g*_C [m^2/s]','particle conservation','electron continuity', 'g*_A [m^2/s]', 'E_{w,A} [V/m]','momentum balance','ExitFlag'};
varTypes2=repmat({'double'},1,length(ErrorNames));
sz2= [1 length(ErrorNames)];
Error = table('Size', sz2, 'VariableTypes',varTypes2, 'VariableNames', ErrorNames);

%%
ctr = 0;
ctr2 = 0;
strt = 1;
% 
% x0=10;
% y0=300;
% width=670;
% height=490;
% set(gcf,'position',[x0,y0,width,height])

exfl = zeros(n_iter, T_iter);

itermat = zeros(n_iter, T_iter);
funcmat = zeros(n_iter, T_iter);

optmat = zeros(n_iter, T_iter);




F_SEE=1; %Include SEE 
F_TEE=1; %Include TEE
F_FEE=0; %Include FEE 
% for F_TEE = 0:1
%     F_SEE = 1-F_TEE;
for te = 1%:length(r_T_wkc)
    for j = 2 %:n_iter
        
        Te = r_Te(2); % Te= 2.5*e; %
        T_wka = r_T_wka(3);
        T_wkc = 3100;% r_T_wkc(te);
        nb = 1*10^22;%r_nb(j);       % electron bulk density 10^22; %
        %n_n = (1-a_iz)/a_iz*ne_0;
        phi_A = 1000; %r_phi_A(j);
        %By = 0;%.7;%r_B(j);
        h = 0.05;%r_h(te);
        C_guess = 5%217.81; %r_C_guess(1);
        
        %
        
        % Flags
        if F_SEE
        SEE_inc = 'yes';
        else
            SEE_inc = 'a posteriori';
        end
 
        if F_TEE
        TEE_inc = 'yes';
        else
        TEE_inc = 'a posteriori';
        end

        if F_FEE
            FEE_inc = 'yes';
        else
            FEE_inc = 'a posteriori';
        end
        
        A_G1 = A_G*F_TEE;%      % Material Constant for shottkey equation
        
        
        
        plasma_properties = {Te, nb,a_iz,Z, m_i};
        design_parameters = {T_wka, T_wkc, E_i, A_G1, h, L, W, E_F};
        
        %% First iteration, forcing E_w = 0
        fix = 0;
        [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By, x,fx,exitflag,initial_state,output] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C, u_ze, d_J, C_guess,  F_SEE, F_TEE, F_FEE,fix)

        output;
        
            % results including non resolved 
            if isreal(x)
            ctr2 = ctr2 + 1;
            InpuValues = {nb, Te/e,phi_A, T_wkc, T_wka, h, d_J, SEE_inc, TEE_inc, FEE_inc};
            Input(ctr2,:) = InpuValues;         
            Output(ctr2,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By};
            Potential(ctr2,:) = {V_C, V_A, E_wc, E_wa, phi_B, phi_D, (phi_B - phi_D) / h};
            IniTable(ctr2,:) = array2table([initial_state phi_A-initial_state(2)*Te/e phi_C-initial_state(1)*Te/e]);
            
            geC_bolz = ge_bolz(nb/2, Te, -V_C)*e;
            SEE_C = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
            TEE_C = schottky(T_wkc, W, E_wc, A_G)*e;
            J = -e*nb*uxe;
            geA_bolz = ge_bolz(nb/2, Te, -V_A)*e;
            SEE_A = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
            TEE_A = schottky(T_wka, W, E_wc, A_G)*e;
            
            Emissions(ctr2,:) = {geC_bolz, geC_em*e, SEE_C, TEE_C, geA_bolz, geA_em*e, SEE_A, TEE_A, J};
            Error(ctr2,:) = array2table([fx exitflag]); %num2cell(fx);
            end

        
        %         
%         E_wc = (E_wc+8.3*10^7)/2
%         geC_em2 = schottky(T_wkc, W, E_wc, A_G1,0);
%         initial_state = [V_C*e/Te, V_A*e/Te,geC_em2, geA_em, E_wc, E_wa, uxe-(geC_em2-geC_em)/(2*nb)];% (E_wc+8.3*10^7)/2
%         
%         [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By, x,fx,exitflag,initial_state,output] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C, u_ze, d_J, C_guess,  F_SEE, F_TEE, F_FEE,0,initial_state);
%       
        %% Second iteration forcing E_w (in) = E_w(output) of previous iteration 
        E_wc1 = E_wc
        geC_em2 = F_TEE*schottky(T_wkc, W, E_wc , A_G1) + F_SEE*SEE(nb/2/Z*sqrt(Te/m_i), E_i, W); %geC_em+schottky(T_wkc, W, E_wc1 , A_G1,0)-schottky(T_wkc, W, E_wc , A_G1,0);
%         
%         E_guess2 = wall_e_field(T_wkc, V_C*e/Te,geC_em2,nb/2, Te)
%         gem_guess2 = F_TEE*schottky(T_wkc, W, E_guess2 , A_G1) + F_SEE*SEE(nb/2/Z*sqrt(Te/m_i), E_i, W);
% 
%         E_guess3 = wall_e_field(T_wkc, V_C*e/Te,gem_guess2 ,nb/2, Te)
%         gem_guess3 = F_TEE*schottky(T_wkc, W, E_guess3, A_G1) + F_SEE*SEE(nb/2/Z*sqrt(Te/m_i), E_i, W);
% 
% lmct = 0;
% 
% while abs(E_guess3-E_guess2)/E_guess2>10^-3
%     
%     if lmct>20
%         disp("Intial guess for E_w did not converge after 20 iterations")
%         break
%     else
%        lmct = lmct+1; 
%     end
%     E_guess2 = wall_e_field(T_wkc, V_C*e/Te,gem_guess3,nb/2, Te)
%     gem_guess2 = F_TEE*schottky(T_wkc, W, E_guess2 , A_G1) + F_SEE*SEE(nb/2/Z*sqrt(Te/m_i), E_i, W);
%     
%     E_guess3 = wall_e_field(T_wkc, V_C*e/Te,gem_guess2 ,nb/2, Te)
%     gem_guess3 = F_TEE*schottky(T_wkc, W, E_guess3, A_G1) + F_SEE*SEE(nb/2/Z*sqrt(Te/m_i), E_i, W);
% 
% end
% 
%     E_wc1 = wall_e_field(T_wkc, V_C*e/Te,gem_guess3 ,nb/2, Te); % E_guess1;%
%     geC_em2 = F_TEE*schottky(T_wkc, W, E_wc1, A_G1) + F_SEE*SEE(nb/2/Z*sqrt(Te/m_i), E_i, W);
% 
%     
        initial_state = [V_C*e/Te, V_A*e/Te,geC_em2, geA_em, E_wc1 , E_wa, uxe-(geC_em2-geC_em)/(2*nb)];% (E_wc+8.3*10^7)/2
        
        [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By, x,fx,exitflag,initial_state,output] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C, u_ze, d_J, C_guess,  F_SEE, F_TEE, F_FEE,E_wc1,initial_state);

        output
            % results including non resolved 
            if isreal(x)
            ctr2 = ctr2 + 1;
            InpuValues = {nb, Te/e,phi_A, T_wkc, T_wka, h, d_J, SEE_inc, TEE_inc, FEE_inc};
            Input(ctr2,:) = InpuValues;         
            Output(ctr2,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By};
            Potential(ctr2,:) = {V_C, V_A, E_wc, E_wa, phi_B, phi_D, (phi_B - phi_D) / h};
            IniTable(ctr2,:) = array2table([initial_state phi_A-initial_state(2)*Te/e phi_C-initial_state(1)*Te/e]);
            
            geC_bolz = ge_bolz(nb/2, Te, -V_C)*e;
            SEE_C = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
            TEE_C = schottky(T_wkc, W, E_wc, A_G)*e;
            J = -e*nb*uxe;
            geA_bolz = ge_bolz(nb/2, Te, -V_A)*e;
            SEE_A = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
            TEE_A = schottky(T_wka, W, E_wc, A_G)*e;
            
            Emissions(ctr2,:) = {geC_bolz, geC_em*e, SEE_C, TEE_C, geA_bolz, geA_em*e, SEE_A, TEE_A, J};
            Error(ctr2,:) = array2table([fx exitflag]); %num2cell(fx);
            end

        %if exitflag<=0
        %% Iterating between E_W is not forced and E_w is forced
        
        now = 0
        for i = 1:OutLim 
            disp("difference with previous E_w: ")
           disp(E_wc1 - E_wc) 
           
           if (fix == -1)&&(exitflag>0)% (abs(E_wc1 - E_wc)/E_wc <epres)&&
               disp("Outer iteration loop converged")
               break            
           end
           disp("fx(1): ")
           disp(fx(1))
        if abs(E_wc1 - E_wc)/E_wc >10*epres%&&(now==0)%&&(exitflag<0)
            fix = E_wc%(E_wc+E_wc1)/2
            E_wc1 = E_wc%(E_wc+E_wc1)/2;
%             now = 0;
%         elseif exist('f2','var')&&(fix==-1)&&now% secant
%             fix = -1  
%             
%             f1 = fx(1)%j_Ewc;
%             x1 = E_wc;
%             
%             delta = f1*(x1-x2)/(f1-f2)
%             
%             E_wc1 = E_wc - delta
%             
%             f2 = f1;
%             x2 = x1;
%         elseif (fix==-1)&(now==0) %if (abs(fx(1)) < epres*E_wc*10^6)&&(fix==-1) % raphson 
%             fix = -1
%             E_wc
%             j_Ewc
%             f2 = fx(1)
%             E_wc1 = E_wc + fx(1)
%             %j_Ewc
%             x2 = E_wc1;
%             now = 1
        else
            fix = -1  
            E_wc1 = E_wc% - fx(1);
        end
            
            geC_em2 = F_TEE*schottky(T_wkc, W, E_wc , A_G1) + F_SEE*SEE(nb/2/Z*sqrt(Te/m_i), E_i, W); %geC_em2 = geC_em+schottky(T_wkc, W, E_wc1 , A_G1)-schottky(T_wkc, W, E_wc , A_G1);
        initial_state = [V_C*e/Te, V_A*e/Te,geC_em2, geA_em, E_wc1 , E_wa, uxe-(geC_em2-geC_em)/(2*nb)];% (E_wc+8.3*10^7)/2
        
        [V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By, x,fx,exitflag,initial_state,output,j_Ewc] = transversal_V3(plasma_properties, design_parameters, phi_A,phi_C, u_ze, d_J, C_guess,  F_SEE, F_TEE, F_FEE,fix,initial_state);
        
        if (i == OutLim)&&((fix > -1)||(exitflag<=0))
            fprintf("Outer iteration loop reached iteration limit (%d), no convergence \n",OutLim)
        end
            
        
        output
            % results including non resolved 
            if isreal(x)
            ctr2 = ctr2 + 1;
            InpuValues = {nb, Te/e,phi_A, T_wkc, T_wka, h, d_J, SEE_inc, TEE_inc, FEE_inc};
            Input(ctr2,:) = InpuValues;         
            Output(ctr2,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By};
            Potential(ctr2,:) = {V_C, V_A, E_wc, E_wa, phi_B, phi_D, (phi_B - phi_D) / h};
            IniTable(ctr2,:) = array2table([initial_state phi_A-initial_state(2)*Te/e phi_C-initial_state(1)*Te/e]);
            
            geC_bolz = ge_bolz(nb/2, Te, -V_C)*e;
            SEE_C = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
            TEE_C = schottky(T_wkc, W, E_wc, A_G)*e;
            J = -e*nb*uxe;
            geA_bolz = ge_bolz(nb/2, Te, -V_A)*e;
            SEE_A = SEE(nb/2/Z*sqrt(Te/m_i), E_i, W)*e;
            TEE_A = schottky(T_wka, W, E_wc, A_G)*e;
            
            Emissions(ctr2,:) = {geC_bolz, geC_em*e, SEE_C, TEE_C, geA_bolz, geA_em*e, SEE_A, TEE_A, J};
            Error(ctr2,:) = array2table([fx exitflag]); %num2cell(fx);
            end

        
        end
        if exitflag>0&&isreal(x)
            exitflag;
            ctr = ctr+1;
            InpuTable(ctr,:) = {nb, Te/e,phi_A, T_wkc, T_wka, h, d_J};
            
            OutpuTable(ctr,:) = {V_C, V_A, geC_em, geA_em, E_wc, E_wa, uxe, phi_B, phi_D, By};
            
            ExiTable(ctr,:) = array2table([exitflag fx]);
        end
%         
%         itermat(j,te) = output.iterations;
%         funcmat(j,te) = output.funcCount;
%         optmat(j,te) = output.firstorderopt;
%         
%         exfl(j,te) = exitflag;
    end
    strt = ctr+1;
end
%end

Input(ctr2+1:end,:)=[]
Output(ctr2+1:end,:)=[]
Potential(ctr2+1:end,:)=[]
Emissions(ctr2+1:end,:)=[]
Error(ctr2+1:end,:)=[]


InpuTable(ctr+1:end,:)=[];
IniTable(ctr2+1:end,:)=[]
OutpuTable(ctr+1:end,:)=[];
ExiTable(ctr+1:end,:)=[];
%% Plot outputs
% for i = 1:7
% figure(i)
% plot(Input{:,4},Output{:,i})
% %plot(Input{:,4},Error{:,i})
% %hold on
% xlabel(InputVarNames(4))
% %ylabel(ErrorNames(i))
% ylabel(OutputVarNames(i))
% 
% end



%% Make Matrices
% 
% %matx = [r_nb,r_Te/e,InpuTable.('Electron density'),InpuTable.('Electrode temperature [eV]')]; 
% label_1 = 'Electron density [m^{-3}]'; % 'Anode potential';%'Magnetic Field'; %
% label_2 = 'Electrode temperature [eV]'; % 'Electrode distance';%
% 
% inmat1 = InpuTable.(label_1); 
% inmat2 = InpuTable.(label_2);
% invec1 = unique(inmat1);
% invec2 = unique(inmat2);
% 
% NT = OrganizeFinds(inmat1,inmat2, InpuTable,OutpuTable)
% %Er = OrganizeFinds(inmat1,inmat2, InpuTable,ExiTable);

% %% Plot with invec2 on the x-axis and color coding representing invec1
% 
% % create color scheme
% cc=jet(n_iter);
% 
% figure(1)
% 
% 
% set(gcf,'position',[100,20,900,640])
% %plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),-OutpuTable.('electron bulk velocity')(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))
%     subplot(2,2,1);
%     hbv = plot(invec2,-NT.Bulkvel');
%     set(hbv, {'color'},num2cell(cc,2))
%     title('Electron Bulk Velocity')
%     xlabel(label_2)
%     ylabel('U_xe [m/s]')
% 
% %figure(2)
% %plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),-OutpuTable.('electron bulk velocity')(strt:ctr).*nb*e,style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))
%     subplot(2,2,2);
%     hcu = plot(invec2,(-NT.Bulkvel.*r_nb'*e)');
%     set(hcu, {'color'},num2cell(cc,2))
%     title('Electron Current')
%     xlabel(label_2)
%     ylabel('I [A/m^2]')
% 
% %figure(3)
% %plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.("Cathode potential")(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))
%     subplot(2,2,3);
%     hcat = plot(invec2,NT.CatPot');
%     set(hcat, {'color'},num2cell(cc,2))
%     title('Cathode Potential Drop')
%     xlabel(label_2)
%     ylabel('Potential Difference [V]')
% 
% %figure(4)
%     subplot(2,2,4);
%     hres = plot(invec2,NT.OuBy');
%     set(hres, {'color'},num2cell(cc,2))
%     %title('Bulk Potentail Drop')
%     title('Magnetic Field')
%     xlabel(label_2)
%     ylabel('B_y [T]')
%     %ylabel('Potential Difference [V]')
%     
% %plot(InpuTable.('Electrode temperature [eV]')(strt:ctr),OutpuTable.("Anode potential")(strt:ctr)-OutpuTable.("Cathode potential")(strt:ctr),style3,'DisplayName', sprintf('n = %.1e m^{-3}', nb),'color',cc(j,:))
% 
% 
% 
% % set ticks labaels for for invec1
% if LOGDENS % in the case invec1= density and logarithmically plotted
%     tiLabls1 = cellfun(@(c) sprintf('%0.1e',c),num2cell(logspace(log10(invec1(1)),log10(invec1(end)),11)),'UniformOutput',false);
% else
%     tiLabls1 = cellfun(@(c) sprintf('%0.1e',c),num2cell(linspace(invec1(1),invec1(end),11)),'UniformOutput',false);
% end
% 
% 
% for i = 1:4
%     subplot(2, 2, mod(i-1, 4)+1)
%     %figure(i)
%     colormap(cc)
%     hc= colorbar;
%     set(hc,'Ticklabels',tiLabls1,'limit',[0 1],'Ticks',linspace(0,1,11))
%     set(hc.Label,'String',label_1)
% end
% 
% 
% %% Plot with invec1 on the x-axis and color coding representing invec2
% 
% 
% %h =  findobj('type','figure');
% %S = 4%length(h);
% figure(2)
% set(gcf,'position',[100,50,900,640])
% 
% c2=jet(length(unique(inmat2)));
%     %figure(S+1)
%     subplot(2,2,1);
%     hrm = plot(invec1,NT.OuBy);
%     set(hrm, {'color'},num2cell(c2,2))
%     xlabel(label_1)
%     %title('Bulk Potentail Drop')
%     title('Magnetic Field')
%     ylabel('B_y [T]')
%     %ylabel('Potential Difference [V]')
% 
%     
%     %figure(S+2)
%     subplot(2,2,2);
%     hcp = plot(invec1,NT.CatPot);
%     set(hcp, {'color'},num2cell(c2,2))
%     title('Cathode Potential Drop')
% xlabel(label_1)
% ylabel('Potentail Drop [V]')
% 
%     subplot(2,2,3);
%     %figure(S+3)
%     hce = plot(invec1,-NT.Bulkvel);
%     set(hce, {'color'},num2cell(c2,2))
% title('Electron Bulk Velocity')
% xlabel(label_1)
% ylabel('U_xe [m/s]')
% 
%     %figure(S+4)
%     subplot(2,2,4);
%     hce = plot(invec1,-NT.Bulkvel.*r_nb'*e);
%     set(hce, {'color'},num2cell(c2,2))
%     title('Electron Current')
% xlabel(label_1)
% ylabel('I [A/m^2]')
% 
% 
% %figure(S+1)
% 
% %figure(S+2)
% % 
% %figure(S+3)
% %title('Cathode Emissions')
% %ylabel('Electron flux [m^{-3}/s]')
% 
% %figure(S+4)
% 
% 
% 
% % Coroize temperaturs
% 
% %Temps = cellfun(@(c) sprintf('%0.1e',c),num2cell([r_Te(1)/e r_Te(T_iter/10:T_iter/10:T_iter)/e]),'UniformOutput',false); % make robost to change
% 
% tiLabls2 = cellfun(@(c) sprintf('%0.1e',c),num2cell(linspace(invec2(1),invec2(end),11)),'UniformOutput',false); % 
% 
% for i = 1:4%S+1:S+4
%     subplot(2, 2, mod(i-1, 4)+1)
%     %figure(i)
%     colormap(c2)
%     hc = colorbar;
%     set(hc,'Ticklabels',tiLabls2,'limit',[0 1],'Ticks',linspace(0,1,11))
%     set(hc.Label,'String',label_2)
%     
% end
% 
% 
% %% 3D plot to check
% S = 2;
% [X,Y] = meshgrid(invec1, invec2);
% C = X.*(Y*e/m_i).^0.5;
% 
%     figure(S+1)
%     set(gcf,'position',[100,20,900,640])
%     subplot(2,2,1);
%     surf(X,Y,NT.OuBy');
%     %view(0,90);
%     colorbar
% 
%     %figure(S+2)
%     subplot(2,2,2);
%     surf(X,Y,NT.CatPot');
%     %view(0,90);
%     colorbar
% 
%     %figure(S+3)
%     subplot(2,2,3);
%     surf(X,Y,-NT.Bulkvel');
%     %view(0,90);
%     colorbar
% 
%     %figure(S+4)
%     subplot(2,2,4);
%     surf(X,Y,(-NT.Bulkvel.*r_nb'*e)');
%     %view(0,90);
%     colorbar
% 
% 
% 
% %figure(S+1)
% subplot(2,2,1);
% title('Magnetic Field')
%     xlabel(label_1)
%     ylabel(label_2)
% 
% %figure(S+2)
% subplot(2,2,2);
% title('Cathode Potential Drop')
%     xlabel(label_1)
%     ylabel(label_2) 
%     
% %figure(S+3)
% subplot(2,2,3);
% title('Electron Bulk Velocity')
%     xlabel(label_1)
%     ylabel(label_2)
% 
% %figure(S+4)
% subplot(2,2,4);
% title('Electron Current')
%     xlabel(label_1)
%     ylabel(label_2)

%% Reorganising outputs into a structure

%load chirp.mat
%sound(y)

function OutStruct = OrganizeFinds(inmat1,inmat2, InpuTable,OutpuTable)
    OutStruct.OuTe = Tab2Mat(inmat1, inmat2, InpuTable.('Electrode temperature [eV]'));
    OutStruct.InNe = Tab2Mat(inmat1, inmat2, InpuTable.('Electron density [m^{-3}]'));
    OutStruct.OuBy = Tab2Mat(inmat1, inmat2, OutpuTable.('Magnetic Field'));
    OutStruct.InH = Tab2Mat(inmat1, inmat2, InpuTable.('Electrode distance'));
    OutStruct.InVA = Tab2Mat(inmat1, inmat2, InpuTable.('Anode potential'));
    OutStruct.CatDrop = Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode sheath potential drop'));
    OutStruct.AnoDrop = Tab2Mat(inmat1, inmat2, OutpuTable.('Anode sheath potential drop'));
    OutStruct.CatEm = Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode emissions'));
    OutStruct.AnoEm = Tab2Mat(inmat1, inmat2, OutpuTable.('Anode emissions'));
    OutStruct.CatEfield = Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode E-field'));
    OutStruct.AnoEfield = Tab2Mat(inmat1, inmat2, OutpuTable.('Anode E-field'));
    OutStruct.Bulkvel =  Tab2Mat(inmat1, inmat2, OutpuTable.('electron bulk velocity'));
    if any(strcmp(OutpuTable.Properties.VariableNames,'Cathode potential'))
        OutStruct.CatPot  =  Tab2Mat(inmat1, inmat2, OutpuTable.('Cathode potential'));
        OutStruct.AnoPot  =  Tab2Mat(inmat1, inmat2, OutpuTable.('Anode potential'));
    %OutStruct.MagPot  = OutStruct.OuBy*u_ze./OutStruct.OuH;
        OutStruct.ResMat  =  OutStruct.AnoPot-OutStruct.CatPot;
    end 
    
    % ResPot = m_e/e*collRate_ei(NT.OuNe/Z,Z,NT.OuTe*e).*NT.Bulkvel.*NT.OuH
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
