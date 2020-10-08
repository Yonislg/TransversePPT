% SChottkey test
% https://books.google.nl/books?id=y0FF19lud5YC&pg=PA6&lpg=PA6&dq=&redir_esc=y#v=onepage&q&f=false
import transversemodel.subfunctions.*;
e = 1.602176634e-19;
a=[1 2 3 4 5]*10^4;
E_w = ([1 2 3 4 5]*10^4).^2
A_G = 1.2016e+06

y1=schottky(1800, 2.5, E_w, A_G)*e
y2=schottky(1800, 3, E_w, A_G)*e

semilogy(a,y1,a,y2)
title('Schottkey Verification using ZrO/W Cathode example')
ylabel('Current density [A/m^2]')
xlabel('(E_W [V/m])^{1/2}') 