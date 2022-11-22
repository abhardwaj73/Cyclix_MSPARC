close all
clear
clc

% Inputs
a0 = 2.67299792458;
n = 13;
m = 13;

%alpha_ext = 0.05;


% 2D hexagonal sheet
%-------------------
a = sqrt(3) * a0;
a1 = [sqrt(3) 1]*a/2;
a2 = [sqrt(3) -1]*a/2;

% FBZ for 2D sheet
%-----------------
K1_hex = [1 sqrt(3)]*2*pi/(sqrt(3)*a);
K2_hex = [1 -sqrt(3)]*2*pi/(sqrt(3)*a);
nd_red_hex = [2/3 1/3;1/3 -1/3; -1/3 -2/3;-2/3 -1/3;-1/3 1/3;1/3 2/3]; % special points in BZ where Dirac cones meet in reduced coordinates
%nd_red_hex = [0 0;1/3 -1/3; 1/2 0];
nd_hex = nd_red_hex(:,1) * K1_hex + nd_red_hex(:,2) * K2_hex;
nd_hex_closed = [nd_hex;nd_hex(1,:)];
qcart_hex = nd_hex;

figure1 = figure;
hold on
box on
scatter(qcart_hex(:,1),qcart_hex(:,2),'r','Marker','o');
plot(nd_hex_closed(:,1),nd_hex_closed(:,2),'r','LineWidth',1,'Linestyle','-.');
%labelpoints(qcart_hex(:,1),qcart_hex(:,2),{'1h', '2h', '3h', '4h', '5h', '6h'});

% Periodic unfolded tube
%-----------------------
d = gcd(n,m);
dR = gcd(2*n+m, 2*m+n);
Ch = a * sqrt(n^2 + m^2 + n*m); % Circumference
T = sqrt(3) * Ch/dR;
N_units = 2 * (n^2 + m^2 + n*m)/dR;
Ch_theta = acos((2*n+m)/(2*sqrt(n^2+m^2+n*m)))*180/pi;
Ch_uvec = (n*a1 + m*a2)/Ch;
T_uvec = (((2*m+n)/dR)*a1 - ((2*n+m)/dR)*a2)/T;

% FBZ for periodic unfolded tube
%-------------------------------
K1_rec = (2*pi/Ch) * Ch_uvec;
K2_rec = (2*pi/T) * T_uvec;
nd_rec = [K1_rec + K2_rec; -K1_rec + K2_rec; -K1_rec - K2_rec; K1_rec - K2_rec]/2;
nd_rec_closed = [nd_rec;nd_rec(1,:)];
L1_r = norm(K1_rec); L2_r = norm(K2_rec);

%R_rec = [sqrt(3)*(n+m)/Ch (m-n)/Ch;(n-m)/Ch 3*(n+m)/T/dR]*a/2;
R_rec = [K1_rec(1)/L1_r K2_rec(1)/L2_r;K1_rec(2)/L1_r K2_rec(2)/L2_r];
q_rec = (R_rec\qcart_hex')';
q_rec(:,1) = q_rec(:,1) + 0.5*L1_r;
q_rec(:,1) = mod(q_rec(:,1),L1_r) - 0.5*L1_r;
q_rec(:,2) = q_rec(:,2) + 0.5*L2_r;
q_rec(:,2) = mod(q_rec(:,2),L2_r) - 0.5*L2_r;
qcart_rec = (R_rec*(q_rec'))'; % convert into cartesian

% To check whether the mapping of unit cell throughout the rec. space puts special points on top of each other
%-------------------------------------------------------------------------------------------------------------
% sz = size(qcart_rec,1);
% Q_rec = zeros(sz*39*39,2);
% count = 1;
% for i = -19:19
% 	for j = -19:19
% 		Q_rec(count:count+5,:) = q_rec + [i*L1_r j*L2_r];
% 		count = count + 6;
% 	end
% end
% Qcart_rec = (R_rec*(Q_rec'))'; 

scatter(qcart_rec(:,1),qcart_rec(:,2),'g','Marker','*');
%scatter(Qcart_rec(:,1),Qcart_rec(:,2),'b','Marker','*');
plot(nd_rec_closed(:,1),nd_rec_closed(:,2),'g','LineWidth',1,'Linestyle','-.');
%labelpoints(qcart_rec(:,1),qcart_rec(:,2),{'1r', '2r', '3r', '4r', '5r', '6r'});

% Parallelogram unit cell with two atoms
%---------------------------------------
flag = 0; 
if m == 0
    M = 1;
    flag = 1;
else
    M_1 = (2*n+m)*d/m/dR;
    lowerb = floor(-m*M_1/N_units + 1); % lower bound
    upperb = ceil(m - m*M_1/N_units - 1); % upper bound
    for q = lowerb:upperb
        M = M_1 + N_units*q/m;
        p = (d + q*n)/m;
        if (round(M) == M && round(p) == p)
            flag = 1;
            break;
        end
    end
end

%M = M + alpha_ext*d*T/2/pi;


% FBZ for the unfolded cyclix (parallelogram) unit cell
%--------------------------------------
K1_par = d*K1_rec - M*K2_rec;
K2_par = (N_units/d)*K2_rec;
nd_par = [K1_par + K2_par; -K1_par + K2_par; -K1_par - K2_par; K1_par - K2_par]/2;
nd_par_closed = [nd_par;nd_par(1,:)];
L1_p = norm(K1_par); L2_p = norm(K2_par);
R_par = [K1_par(1)/L1_p K2_par(1)/L2_p;K1_par(2)/L1_p K2_par(2)/L2_p];
q_par = (R_par\qcart_hex')';
q_par(:,1) = q_par(:,1) + 0.5*L1_p;
q_par(:,1) = mod(q_par(:,1),L1_p) - 0.5*L1_p;
q_par(:,2) = q_par(:,2) + 0.5*L2_p;
q_par(:,2) = mod(q_par(:,2),L2_p) - 0.5*L2_p;
qcart_par = (R_par*(q_par'))'; % convert into cartesian

scatter(qcart_par(:,1),qcart_par(:,2),'b','Marker','d');
plot(nd_par_closed(:,1),nd_par_closed(:,2),'b','LineWidth',1,'Linestyle','-');
%labelpoints(qcart_par(:,1),qcart_par(:,2),{'1p', '2p', '3p', '4p', '5p', '6p'});

% save the plot
%--------------
xlim([-32 32]);
ylim([-32 32]);
xlabel('$\textbf{k}_x$','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\textbf{k}_y$','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'NT_FBZ','epsc')
hold off

% FBZ for the unfolded cyclix tube in helical coordinate system
K1_hel = [1 0]*d*2*pi/Ch;
K2_hel = [0 1]*(N_units/d)*2*pi/T;
nd_hel = [K1_hel + K2_hel; -K1_hel + K2_hel; -K1_hel - K2_hel; K1_hel - K2_hel]/2;
nd_hel_closed = [nd_hel;nd_hel(1,:)];
L1_hel = norm(K1_hel); L2_hel = norm(K2_hel);
%R_hel = [1 -M*Ch/d/T;0 1];
R_hel = [1 0;M*Ch/d/T 1];
%q_hel = (R_hel*(q_rec'))';
q_hel = (R_hel*(R_rec\(R_par*q_par')))';
q_hel(:,1) = q_hel(:,1) + 0.5*L1_hel;
q_hel(:,1) = mod(q_hel(:,1),L1_hel) - 0.5*L1_hel;
q_hel(:,2) = q_hel(:,2) + 0.5*L2_hel;
q_hel(:,2) = mod(q_hel(:,2),L2_hel) - 0.5*L2_hel;

q_hel./[L1_hel L2_hel]

figure2 = figure;
hold on
box on

scatter(q_hel(:,1),q_hel(:,2),'b','Marker','d');
plot(nd_hel_closed(:,1),nd_hel_closed(:,2),'b','LineWidth',1,'Linestyle','-');
%labelpoints(q_hel(:,1),q_hel(:,2),{'1', '2', '3', '4', '5', '6'});

xlim([-32 32]);
ylim([-32 32]);
xlabel('$\textbf{k}_{\tilde{x}}$','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\textbf{k}_{\tilde{y}}$','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

set(figure2,'Units','Inches');
pos = get(figure2,'Position');
set(figure2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'NT_FBZ_helical','epsc')
hold off