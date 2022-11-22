close all
clear
clc

% Carbon nanotube 1
%----------------

% Radius = 11.80 Bohr, mesh-size=0.125, N_c = 16, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

tw_d1 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr ->rad/nm

Epl_d1 = [-2.2865311400E+01;-2.2865247940E+01;-2.2865058001E+01;-2.2864741106E+01;-2.2864296698E+01] * 16*27.2114 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_d1 = Epl_d1 - Epl_d1(1);

tw_d1 = tw_d1(2:end);
rEpl_d1 = rEpl_d1(2:end);

log_tw_d1 = log10(tw_d1);
log_rEpl_d1 = log10(rEpl_d1);
P_Efit_d1 = polyfit(log_tw_d1,log_rEpl_d1,1);
Efit_d1 = P_Efit_d1(1)*log_tw_d1 + P_Efit_d1(2);
antilog_Efit_d1= 10.^(Efit_d1);


% Carbon nanotube 2
%----------------

% Radius = 14 Bohr, mesh-size=0.125, N_c = 19, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

tw_d2 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr ->rad/nm

Epl_d2 = [-2.2867450306E+01;-2.2867360928E+01;-2.2867092846E+01;-2.2866645245E+01;-2.2866017407E+01] * 27.2114*19 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_d2 = Epl_d2 - Epl_d2(1);

tw_d2 = tw_d2(2:end);
rEpl_d2 = rEpl_d2(2:end);

log_tw_d2 = log10(tw_d2);
log_rEpl_d2 = log10(rEpl_d2);
P_Efit_d2 = polyfit(log_tw_d2,log_rEpl_d2,1);
Efit_d2 = P_Efit_d2(1)*log_tw_d2 + P_Efit_d2(2);
antilog_Efit_d2= 10.^(Efit_d2);

% Carbon nanotube 3
%----------------

% Radius = 16.18 Bohr, mesh-size=0.125, N_c = 22, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

tw_d3 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr ->rad/nm

Epl_d3 = [-2.2868769901E+01;-2.2868650159E+01;-2.2868290674E+01;-2.2867689807E+01;-2.2866846540E+01] * 27.2114*22 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_d3 = Epl_d3 - Epl_d3(1);

tw_d3 = tw_d3(2:end);
rEpl_d3 = rEpl_d3(2:end);

log_tw_d3 = log10(tw_d3);
log_rEpl_d3 = log10(rEpl_d3);
P_Efit_d3 = polyfit(log_tw_d3,log_rEpl_d3,1);
Efit_d3 = P_Efit_d3(1)*log_tw_d3 + P_Efit_d3(2);
antilog_Efit_d3= 10.^(Efit_d3);

% Carbon nanotube 4
%----------------

% Radius = 18.38 Bohr, mesh-size=0.125, N_c = 25, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

tw_d4 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr ->rad/nm

Epl_d4 = [-2.2869641143E+01;-2.2869486655E+01;-2.2869021988E+01;-2.2868245291E+01;-2.2867154569E+01] * 27.2114*25 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_d4 = Epl_d4 - Epl_d4(1);

tw_d4 = tw_d4(2:end);
rEpl_d4 = rEpl_d4(2:end);

log_tw_d4 = log10(tw_d4);
log_rEpl_d4 = log10(rEpl_d4);
P_Efit_d4 = polyfit(log_tw_d4,log_rEpl_d4,1);
Efit_d4 = P_Efit_d4(1)*log_tw_d4 + P_Efit_d4(2);
antilog_Efit_d4= 10.^(Efit_d4);


% Carbon nanotube 5
%----------------

% Radius = 20.57 Bohr, mesh-size=0.125, N_c = 28, N_eta = 15, H = 7.98 Bohr

H = 7.98;
L = 4.607255148133214;

tw_d5 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr ->rad/nm

Epl_d5 = [-2.2870246344E+01;-2.2870052800E+01;-2.2869469904E+01;-2.2868494282E+01;-2.2867124039E+01] * 27.2114*28 /(H/18.8973 * L/18.8973) ; % input is in Ha -> eV/nm^2

rEpl_d5 = Epl_d5 - Epl_d5(1);

tw_d5 = tw_d5(2:end);
rEpl_d5 = rEpl_d5(2:end);

log_tw_d5 = log10(tw_d5);
log_rEpl_d5 = log10(rEpl_d5);
P_Efit_d5 = polyfit(log_tw_d5,log_rEpl_d5,1);
Efit_d5 = P_Efit_d5(1)*log_tw_d5 + P_Efit_d5(2);
antilog_Efit_d5= 10.^(Efit_d5);


% Twist energy Vs twist plot
%------------

figure1 = figure;
box on

hold on
scatter(tw_d1,rEpl_d1,'b','Marker','o')
scatter(tw_d2,rEpl_d2,'r','Marker','*')
scatter(tw_d3,rEpl_d3,'g','Marker','sq')
scatter(tw_d4,rEpl_d4,'c','Marker','d')
scatter(tw_d5,rEpl_d5,'m','Marker','+')
plot(tw_d1,antilog_Efit_d1,'color','b','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d2,antilog_Efit_d2,'color','r','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d3,antilog_Efit_d3,'color','g','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d4,antilog_Efit_d4,'color','c','Markersize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d5,antilog_Efit_d5,'color','m','Markersize',10,'LineWidth',1,'LineStyle','-.')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'xscale','log','yscale','log','XTick',[0.01 0.02 0.04],'XMinorTick','on','YTick',[1e-2 1e-1 1e0],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
%set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
set(gca,'YTickLabel',{'$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
xlim([0.008 0.04]);
ylim([1e-2 1]);

xlabel('$\alpha$ (rad/nm)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\mathcal{E}_{twist}$ (eV/nm$^2$)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

%ylabel('Relative percentage error (%)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
legend1 = legend('$d = 1.25$ nm','$d = 1.48$ nm','$d = 1.71$ nm','$d = 1.94$ nm','$d = 2.18$ nm','Location','Southeast');
set(legend1,'fontsize',14,'FontName','Times New Roman','Interpreter','LaTex');

str = 'Slope = 2.0';
dim = [.2 .58 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'twistenergyVstwist','epsc')
hold off
