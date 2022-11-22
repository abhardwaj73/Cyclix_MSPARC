close all
clc
clear

n = 6;
m = 3;
a = 2.67;

A = sqrt(3)*a;
b1 = [1/sqrt(3) 1]*2*pi/A; b2 = [1/sqrt(3) -1]*2*pi/A;

% High symmetry points in Brillouin zone of hexagonal lattice (vertices)
xx = [1/3 2/3 1/3 -1/3 -2/3 -1/3];
yy = [-1/3 1/3 2/3 1/3 -1/3 -2/3];

X = xx'*b1;
Y = yy'*b2;
Z = X+Y;
ZZ = [Z;Z(1,:)];
figure1 = figure;
hold on
scatter(Z(:,1),Z(:,2),'r','Marker','o');
plot(ZZ(:,1),ZZ(:,2),'r','LineWidth',1,'Linestyle','-.');

% Reciprocal space lattice vectors

d = gcd(n,m);
dR = gcd(2*n+m, 2*m+n);
Ch = A * sqrt(n^2 + m^2 + n*m); % Circumference
T = sqrt(3) * Ch/dR;
N_units = 2 * (n^2 + m^2 + n*m)/dR;

K1 = (d*A*A/2/Ch/Ch) * ( (2*n + m) * b1 + (2*m + n) * b2 );
K2 = (3*N_units*A*A/2/d/T/T/dR) * ( m * b1 - n * b2 );

f1 = @(x) -(K1(1)/K1(2)) * x + (K1(1)^2 + K1(2)^2)/2/K1(2);
f2 = @(x) -(K2(1)/K2(2)) * x + (K2(1)^2 + K2(2)^2)/2/K2(2);
f3 = @(x) -(K1(1)/K1(2)) * x - (K1(1)^2 + K1(2)^2)/2/K1(2);
f4 = @(x) -(K2(1)/K2(2)) * x - (K2(1)^2 + K2(2)^2)/2/K2(2);

fun1 = @(x) f1(x) - f2(x);
fun2 = @(x) f2(x) - f3(x);
fun3 = @(x) f3(x) - f4(x);
fun4 = @(x) f4(x) - f1(x);

xints_1 = fzero(fun1,1); yints_1 = f1(xints_1);
xints_2 = fzero(fun2,1); yints_2 = f2(xints_2);
xints_3 = fzero(fun3,1); yints_3 = f3(xints_3);
xints_4 = fzero(fun4,1); yints_4 = f4(xints_4);
xints = [xints_1 xints_2 xints_3 xints_4 xints_1];
yints = [yints_1 yints_2 yints_3 yints_4 yints_1];


% Projection of high-symmetry points on BZ of (n,m) tube
L1 = norm(K1); L2 = norm(K2);
K1_hat = K1/L1; K2_hat = K2/L2;
PX = sum(Z .* repmat(K1_hat,size(Z,1),1),2);
PY = sum(Z .* repmat(K2_hat,size(Z,1),1),2);
PX = PX + 0.5*L1;
PY = PY + 0.5*L2;
PX = mod(PX,L1) - 0.5*L1;
PY = mod(PY,L2) - 0.5*L2;
PZ = PX*K1_hat + PY*K2_hat;
scatter(PZ(:,1),PZ(:,2),'b','Marker','*');

PX_hat = PX/L1;
PY_hat = PY/L2;
fprintf('Reduced coordinates of high-symmetry points are: \n');
[PX_hat PY_hat]

% Plotting of axes
N1_x = [-K1(1)/2 K1(1)/2]; N1_y = [-K1(2)/2 K1(2)/2];
plot(N1_x,N1_y,'k','LineWidth',1.5,'Linestyle',':');

N2_x = [-K2(1)/2 K2(1)/2]; N2_y = [-K2(2)/2 K2(2)/2];
plot(N2_x,N2_y,'c','LineWidth',1.5,'Linestyle',':');

plot (xints,yints,'b','LineWidth',1,'Linestyle','--');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'BZ_plot','epsc')
hold off
