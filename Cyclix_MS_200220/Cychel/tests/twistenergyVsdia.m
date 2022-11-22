close all
clear
clc

% Carbon nanotube 1
%----------------

% twist = 0.0 rad/Bohr, mesh-size=0.125, N_c = 16, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;
Go = [16;19;22;25;28];

D_tw1 = [1.25;1.48;1.71;1.94;2.18]; % in nm

Epl_tw1 = [-2.2865311400E+01;-2.2867450306E+01;-2.2868769901E+01;-2.2869641143E+01;-2.2870246344E+01] * 27.2114 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2
%Epl_tw1 = Epl_tw1.*(pi*D_tw1./Go);
% rEpl_tw1 = Epl_tw1;
% 
% %D_tw1 = D_tw1(2:end);
% %rEpl_tw1 = rEpl_tw1(2:end);
% 
% log_D_tw1 = log10(D_tw1);
% log_rEpl_tw1 = log10(rEpl_tw1);
% P_Efit_tw1 = polyfit(log_D_tw1,log_rEpl_tw1,1);
% Efit_tw1 = P_Efit_tw1(1)*log_D_tw1 + P_Efit_tw1(2);
% antilog_Efit_tw1= 10.^(Efit_tw1);


% Carbon nanotube 2
%----------------

% twist = 0.0005 rad/Bohr, mesh-size=0.125, N_c = 19, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

D_tw2 = [1.25;1.48;1.71;1.94;2.18]; % in nm

Epl_tw2 = [-2.2865247940E+01;-2.2867360928E+01;-2.2868650159E+01;-2.2869486655E+01;-2.2870052800E+01] * 27.2114 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2
%Epl_tw2 = Epl_tw2.*(pi*D_tw2./Go);

rEpl_tw2 = Epl_tw2 - Epl_tw1;

% D_tw2 = D_tw2(2:end);
% rEpl_tw2 = rEpl_tw2(2:end);

log_D_tw2 = log10(D_tw2);
log_rEpl_tw2 = log10(rEpl_tw2);
P_Efit_tw2 = polyfit(log_D_tw2,log_rEpl_tw2,1);
Efit_tw2 = P_Efit_tw2(1)*log_D_tw2 + P_Efit_tw2(2);
antilog_Efit_tw2= 10.^(Efit_tw2);

% Carbon nanotube 3
%----------------

% twist = 0.0010 rad/Bohr, mesh-size=0.125, N_c = 22, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

D_tw3 = [1.25;1.48;1.71;1.94;2.18]; % in nm

Epl_tw3 = [-2.2865058001E+01;-2.2867092846E+01;-2.2868290674E+01;-2.2869021988E+01;-2.2869469904E+01] * 27.2114 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_tw3 = Epl_tw3 - Epl_tw1;

% D_tw3 = D_tw3(2:end);
% rEpl_tw3 = rEpl_tw3(2:end);

log_D_tw3 = log10(D_tw3);
log_rEpl_tw3 = log10(rEpl_tw3);
P_Efit_tw3 = polyfit(log_D_tw3,log_rEpl_tw3,1);
Efit_tw3 = P_Efit_tw3(1)*log_D_tw3 + P_Efit_tw3(2);
antilog_Efit_tw3= 10.^(Efit_tw3);

% Carbon nanotube 4
%----------------

% twist = 0.0015 rad/Bohr, mesh-size=0.125, N_c = 25, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

D_tw4 = [1.25;1.48;1.71;1.94;2.18]; % in nm

Epl_tw4 = [-2.2864741106E+01;-2.2866645245E+01;-2.2867689807E+01;-2.2868245291E+01;-2.2868494282E+01] * 27.2114 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_tw4 = Epl_tw4 - Epl_tw1;

% D_tw4 = D_tw4(2:end);
% rEpl_tw4 = rEpl_tw4(2:end);

log_D_tw4 = log10(D_tw4);
log_rEpl_tw4 = log10(rEpl_tw4);
P_Efit_tw4 = polyfit(log_D_tw4,log_rEpl_tw4,1);
Efit_tw4 = P_Efit_tw4(1)*log_D_tw4 + P_Efit_tw4(2);
antilog_Efit_tw4= 10.^(Efit_tw4);


% Carbon nanotube 5
%----------------

% twist = 0.0020 rad/Bohr, mesh-size=0.125, N_c = 28, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

D_tw5 = [1.25;1.48;1.71;1.94;2.18]; % in nm

Epl_tw5 = [-2.2864296698E+01;-2.2866017407E+01;-2.2866846540E+01;-2.2867154569E+01;-2.2867124039E+01] * 27.2114 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_tw5 = Epl_tw5 - Epl_tw1;

% D_tw5 = D_tw5(2:end);
% rEpl_tw5 = rEpl_tw5(2:end);

log_D_tw5 = log10(D_tw5);
log_rEpl_tw5 = log10(rEpl_tw5);
P_Efit_tw5 = polyfit(log_D_tw5,log_rEpl_tw5,1);
Efit_tw5 = P_Efit_tw5(1)*log_D_tw5 + P_Efit_tw5(2);
antilog_Efit_tw5= 10.^(Efit_tw5);


% Twist energy Vs twist plot
%------------

figure1 = figure;
box on

hold on
%scatter(D_tw1,rEpl_tw1,'b','Marker','o')
scatter(D_tw2,rEpl_tw2,'r','Marker','*')
scatter(D_tw3,rEpl_tw3,'g','Marker','sq')
scatter(D_tw4,rEpl_tw4,'c','Marker','d')
scatter(D_tw5,rEpl_tw5,'m','Marker','+')
%plot(D_tw1,antilog_Efit_tw1,'color','b','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw2,antilog_Efit_tw2,'color','r','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw3,antilog_Efit_tw3,'color','g','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw4,antilog_Efit_tw4,'color','c','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw5,antilog_Efit_tw5,'color','m','Markersize',10,'LineWidth',1,'LineStyle','-.')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'xscale','log','yscale','log','XTick',[1.25 1.50 1.75 2.00 2.25],'XMinorTick','on','YTick',[1e-2 1e-1 1e0],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
%set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
set(gca,'XTickLabel',{'$1.25$', '$1.50$', '$1.75$', '$2.00$','$2.25$'});
set(gca,'YTickLabel',{'$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
xlim([1.24 2.25]);
ylim([1e-2 1]);

xlabel('$d$ (nm)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\mathcal{E}_{twist}$ (eV/nm$^2$)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

%ylabel('Relative percentage error (%)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
legend1 = legend('$\alpha = 0.009$','$\alpha = 0.019$','$\alpha = 0.028$','$\alpha = 0.037$','Location','Southeast');
set(legend1,'fontsize',12,'FontName','Times New Roman','Interpreter','LaTex');

str = 'Slope = 2.0';
dim = [.2 .58 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'twistenergyVsdia','epsc')
hold off
