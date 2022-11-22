close all
clear
clc

% Carbon nanotube 1
%-------------

% Radius = 11.80 Bohr, mesh-size=0.125, N_c = 16, N_eta = 15

tw_d1 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr -> rad/nm

bg_d1 = [0.018511693532003; 0.018731759239892;0.019377029056173;0.020407489638717;0.021769085816117] * 27.2114; % input is in Ha
rbg_d1 = bg_d1 - bg_d1(1);

tw_d1 = tw_d1(2:end);
rbg_d1 = rbg_d1(2:end);

log_tw_d1 = log10(tw_d1);
log_rbg_d1 = log10(rbg_d1);
P_bgfit_d1 = polyfit(log_tw_d1,log_rbg_d1,1);
bgfit_d1 = P_bgfit_d1(1)*log_tw_d1 + P_bgfit_d1(2);
antilog_bgfit_d1= 10.^(bgfit_d1);

% Carbon nanotube 2
%-------------

% Radius = 14 Bohr, mesh-size=0.125, N_c = 19, N_eta = 15

tw_d2 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr -> rad/nm

bg_d2 = [0.016127761990841;0.016485113606966;0.017513728381320;0.019106225914502;0.021136617000022] * 27.2114; % input is in Ha
rbg_d2 = bg_d2 - bg_d2(1);

tw_d2 = tw_d2(2:end);
rbg_d2 = rbg_d2(2:end);

log_tw_d2 = log10(tw_d2);
log_rbg_d2 = log10(rbg_d2);
P_bgfit_d2 = polyfit(log_tw_d2,log_rbg_d2,1);
bgfit_d2 = P_bgfit_d2(1)*log_tw_d2 + P_bgfit_d2(2);
antilog_bgfit_d2= 10.^(bgfit_d2);

% Carbon nanotube 3
%-------------

% Radius = 16.18 Bohr, mesh-size=0.125, N_c = 22, N_eta = 15

tw_d3 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr -> rad/nm

bg_d3 = [0.014244279484005;0.014784649276933;0.016299189416469;0.018552288466811;0.021312958621136] * 27.2114; % input is in Ha
rbg_d3 = bg_d3 - bg_d3(1);

tw_d3 = tw_d3(2:end);
rbg_d3 = rbg_d3(2:end);

log_tw_d3 = log10(tw_d3);
log_rbg_d3 = log10(rbg_d3);
P_bgfit_d3 = polyfit(log_tw_d3,log_rbg_d3,1);
bgfit_d3 = P_bgfit_d3(1)*log_tw_d3 + P_bgfit_d3(2);
antilog_bgfit_d3= 10.^(bgfit_d3);

% Carbon nanotube 4
%-------------

% Radius = 18.38 Bohr, mesh-size=0.125, N_c = 25, N_eta = 15

tw_d4 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr -> rad/nm

bg_d4 = [0.012735314357575;0.013508510320171;0.015600962935476;0.018574607450225;0.022079404854676] * 27.2114; % input is in Ha
rbg_d4 = bg_d4 - bg_d4(1);

tw_d4 = tw_d4(2:end);
rbg_d4 = rbg_d4(2:end);

log_tw_d4 = log10(tw_d4);
log_rbg_d4 = log10(rbg_d4);
P_bgfit_d4 = polyfit(log_tw_d4,log_rbg_d4,1);
bgfit_d4 = P_bgfit_d4(1)*log_tw_d4 + P_bgfit_d4(2);
antilog_bgfit_d4= 10.^(bgfit_d4);

% Carbon nanotube 5
%-------------

% Radius = 20.57 Bohr, mesh-size=0.125, N_c = 28, N_eta = 15

tw_d5 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr -> rad/nm

bg_d5 = [0.011505752339886;0.012563917483190;0.015307620780065;0.019025146991836;0.023258787736912] * 27.2114; % input is in Ha
rbg_d5 = bg_d5 - bg_d5(1);

tw_d5 = tw_d5(2:end);
rbg_d5 = rbg_d5(2:end);

log_tw_d5 = log10(tw_d5);
log_rbg_d5 = log10(rbg_d5);
P_bgfit_d5 = polyfit(log_tw_d5,log_rbg_d5,1);
bgfit_d5 = P_bgfit_d5(1)*log_tw_d5 + P_bgfit_d5(2);
antilog_bgfit_d5= 10.^(bgfit_d5);
 
% Bandgap Vs twist plot
%------------

figure1 = figure;
box on

hold on
scatter(tw_d1,rbg_d1,'b','Marker','o');
scatter(tw_d2,rbg_d2,'r','Marker','*');
scatter(tw_d3,rbg_d3,'g','Marker','sq');
scatter(tw_d4,rbg_d4,'c','Marker','d');
scatter(tw_d5,rbg_d5,'m','Marker','+');
plot(tw_d1,antilog_bgfit_d1,'b','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d2,antilog_bgfit_d2,'r','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d3,antilog_bgfit_d3,'g','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d4,antilog_bgfit_d4,'c','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(tw_d5,antilog_bgfit_d5,'m','MarkerSize',10,'LineWidth',1,'LineStyle','-.')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'xscale','log','yscale','log','XTick',[0.01 0.02 0.04],'YTick',[0.005 0.01 0.02 0.04 0.08 0.16 0.32],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.3f'))
set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
set(gca,'YTickLabel',{'$0.005$', '$0.010$', '$0.020$', '$0.040$', '$0.080$', '$0.160$', '$0.320$'});
xlim([0.008 0.04]);
ylim([0.005 0.32]);
xlabel('$\alpha$ (rad/nm)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\Delta_{Bandgap}$ (eV)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

%ylabel('Relative percentage error (%)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
legend1 = legend('$d = 1.25$ nm','$d = 1.48$ nm','$d = 1.71$ nm','$d = 1.94$ nm','$d = 2.18$ nm','Location','Southeast');
set(legend1,'fontsize',14,'FontName','Times New Roman','Interpreter','LaTex');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'BandgapVstwist','epsc')
hold off
