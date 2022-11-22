% evaluate stress from energies at different close-by deformations
clc
clear

natm = 2;
go_cyc = 6;
W = [-5.70196269 -5.7019906666 -5.7018627037 -5.7017432555 -5.7014766530]*natm*go_cyc;

%W = [-5.51198481 -5.514187715 -5.52021020 -5.52294025 -5.52511367]*natm*go_cyc;
dW = W - W(3);
%F = [1.510 1.511 1.51374771440196 1.515 1.516];
F = [1.49 1.50 1.51374771440196 1.52 1.53];
dF = (F - F(3))/F(3);

pp=spline(dF,dW);
p_der=fnder(pp,1);
y_prime=ppval(p_der,dF);

stress = y_prime(3)/F(3)