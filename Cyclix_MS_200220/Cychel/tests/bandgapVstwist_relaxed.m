close all
clear
clc

% Carbon nanotube 1
%-------------

% Radius = 11.80 Bohr, mesh-size=0.125, N_c = 16, N_eta = 15

tw_d1 = [0;0.0005;0.0010;0.0015;0.0020] * 18.8973; % input is in rad/Bohr -> rad/nm

bg_d1 = [0.020710458913380;0.019966645387949;0.020509247580111;0.020774931162439;0.021610736102809] * 27.2114; % input is in Ha
rbg_d1 = bg_d1;% - bg_d1(1);

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

bg_d2 = [0.016127761990841;0.016250810623960;0.016854439664971;0.017860074649181;0.019188206311417] * 27.2114; % input is in Ha
rbg_d2 = bg_d2;% - bg_d2(1);

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

bg_d3 = [0.014244279484005;0.014428264016800;0.015365800507213;0.016805542603472;0.018669137603032] * 27.2114; % input is in Ha
rbg_d3 = bg_d3;% - bg_d3(1);

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

bg_d4 = [0.012735314357575;0.013293357357443;0.014362331090117;0.016315663357158;0.018747434493366] * 27.2114; % input is in Ha




rbg_d4 = bg_d4;% - bg_d4(1);

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

bg_d5 = [0.011505752339886;0.012227602187744;0.013775022316560;0.016251088786197;0.019778602231705] * 27.2114; % input is in Ha
rbg_d5 = bg_d5;% - bg_d5(1);

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
set(gca,'xscale','log','yscale','log','XTick',[0.01 0.02 0.04],'YTick',[0.3 0.35 0.40 0.45 0.50 0.55 0.60 0.65],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.2f'))
set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
set(gca,'YTickLabel',{'$0.30$', '$0.35$', '$0.40$', '$0.45$', '$0.50$', '$0.55$', '$0.60$', '$0.65$'});
xlim([0.008 0.04]);
ylim([0.30 0.65]);
xlabel('$\alpha$ (rad/nm)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$E_G$ (eV)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

%ylabel('Relative percentage error (%)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
legend1 = legend('$d = 1.25$ nm','$d = 1.48$ nm','$d = 1.71$ nm','$d = 1.94$ nm','$d = 2.18$ nm','Location','Southeast');
set(legend1,'fontsize',14,'FontName','Times New Roman','Interpreter','LaTex');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'BandgapVstwist_relaxed','epsc')
hold off
