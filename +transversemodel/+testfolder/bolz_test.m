clear all
close all
global m_e e;
import transversemodel.subfunctions.*;

m_e = 9.1093837015e-31;
e = 1.602176634e-19;

%Potential drop in [V]
x = linspace(0,7);%0 1 2 3 4 5 6 7];


y1 = zeros(3,100);
y2=y1;
% for i = 1:3
%     for j = 1:3
%         y(3*(i-1)+j,:)= ge_bolz((2*j-1)*10^22,i*e,x);
%     end
% end
figure(1)
title("Characteristic Curve for T_e = 1 eV")%Different n_{e0}")
ylabel("Particle Flux [m^{-3}s^{-1}]")
xlabel("Potential Drop [V]") 
hold on
for j = 1:3
     y1(j,:)= ge_bolz((2*j-1)*10^22,1*e,x);
     plot(x,y1(j,:),'DisplayName',strcat("n_{e0} = ",int2str(2*j-1),"\cdot10^{22}m^{-3}"))
end
legend
H = findobj('type','legend')


%figure(1)
%plot(x,y1)
%legend("1","3","5")
figure(2)
title("Characteristic Curve for n_{e0} = 10^{22} m^{-3}")%Different T_e")
ylabel("Particle Flux [m^{-3}s^{-1}]")
xlabel("Potential Drop [V]") 
hold on
for i = 1:3
     y2(i,:)= ge_bolz(10^22, i*e,x);
     plot(x,y2(i,:),'DisplayName',strcat("T_e = ",int2str(i)," eV"))
end
legend
H2 = findobj(figure(2),'type','legend')



%ax = findobj(figure(1),'Type','Axes');
figure(3)
plot(conv(x, [0.5 0.5], 'valid'), diff(y1,1,2))
title("Numerical Derivative for T_e = 1 eV")% various n_{e0}")
ylabel("\Delta Flux / \Delta V")% with \Delta V = 1")
xlabel("Potential Drop [V]") 
legend(H.String)

figure(4)
plot(conv(x, [0.5 0.5], 'valid'), diff(y2,1,2))
title("Numerical Derivative for  n_{e0} = 10^{22} m^{-3}")%various T_e")
ylabel("\Delta Flux / \Delta V")% with \Delta V = 1")
xlabel("Potential Drop [V]") 
legend(H2.String)