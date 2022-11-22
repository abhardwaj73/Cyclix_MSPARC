close all
clear
clc

% Carbon nanotube 1
%-------------

% twist = 0.0 rad/Bohr, mesh-size=0.125, N_c = 16, N_eta = 15

D_tw1 = [1.25;1.48;1.71;1.94;2.18]; % in nm

bg_tw1 = [0.018511693532003;0.016127761990841;0.014244279484005;0.012735314357575;0.011505752339886] * 27.2114; % input is in Ha
% rbg_tw1 = bg_tw1 - bg_tw1(1);
% 
% %D_tw1 = D_tw1(2:end);
% %rbg_tw1 = rbg_tw1(2:end);
% 
% log_D_tw1 = log10(D_tw1);
% log_rbg_tw1 = log10(rbg_tw1);
% P_bgfit_tw1 = polyfit(log_D_tw1,log_rbg_tw1,1);
% bgfit_tw1 = P_bgfit_tw1(1)*log_D_tw1 + P_bgfit_tw1(2);
% antilog_bgfit_tw1= 10.^(bgfit_tw1);

% Carbon nanotube 2
%-------------

% twist = 0.0005 rad/Bohr, mesh-size=0.125, N_c = 19, N_eta = 15

D_tw2 = [1.25;1.48;1.71;1.94;2.18]; % in nm

bg_tw2 = [0.018731759239892;0.016485113606966;0.014784649276933;0.013508510320171;0.012563917483190] * 27.2114; % input is in Ha
rbg_tw2 = bg_tw2 - bg_tw1;

%D_tw2 = D_tw2(2:end);
%rbg_tw2 = rbg_tw2(2:end);

log_D_tw2 = log10(D_tw2);
log_rbg_tw2 = log10(rbg_tw2);
P_bgfit_tw2 = polyfit(log_D_tw2,log_rbg_tw2,1);
bgfit_tw2 = P_bgfit_tw2(1)*log_D_tw2 + P_bgfit_tw2(2);
antilog_bgfit_tw2= 10.^(bgfit_tw2);

% Carbon nanotube 3
%-------------

% twist = 0.0010 rad/Bohr, mesh-size=0.125, N_c = 22, N_eta = 15

D_tw3 = [1.25;1.48;1.71;1.94;2.18]; % in nm

bg_tw3 = [0.019377029056173;0.017513728381320;0.016299189416469;0.015600962935476;0.015307620780065] * 27.2114; % input is in Ha
rbg_tw3 = bg_tw3 - bg_tw1;

%D_tw3 = D_tw3(2:end);
%rbg_tw3 = rbg_tw3(2:end);

log_D_tw3 = log10(D_tw3);
log_rbg_tw3 = log10(rbg_tw3);
P_bgfit_tw3 = polyfit(log_D_tw3,log_rbg_tw3,1);
bgfit_tw3 = P_bgfit_tw3(1)*log_D_tw3 + P_bgfit_tw3(2);
antilog_bgfit_tw3= 10.^(bgfit_tw3);

% Carbon nanotube 4
%-------------

% twist = 0.0015 rad/Bohr, mesh-size=0.125, N_c = 25, N_eta = 15

D_tw4 = [1.25;1.48;1.71;1.94;2.18]; % in nm

bg_tw4 = [0.020407489638717;0.019106225914502;0.018552288466811;0.018574607450225;0.019025146991836] * 27.2114; % input is in Ha
rbg_tw4 = bg_tw4 - bg_tw1;

%D_tw4 = D_tw4(2:end);
%rbg_tw4 = rbg_tw4(2:end);

log_D_tw4 = log10(D_tw4);
log_rbg_tw4 = log10(rbg_tw4);
P_bgfit_tw4 = polyfit(log_D_tw4,log_rbg_tw4,1);
bgfit_tw4 = P_bgfit_tw4(1)*log_D_tw4 + P_bgfit_tw4(2);
antilog_bgfit_tw4= 10.^(bgfit_tw4);

% Carbon nanotube 5
%-------------

% twist = 0.0020 rad/Bohr, mesh-size=0.125, N_c = 28, N_eta = 15

D_tw5 = [1.25;1.48;1.71;1.94;2.18]; % in nm

bg_tw5 = [0.021769085816117;0.021136617000022;0.021312958621136;0.022079404854676;0.023258787736912] * 27.2114; % input is in Ha
rbg_tw5 = bg_tw5- bg_tw1;

%D_tw5 = D_tw5(2:end);
%rbg_tw5 = rbg_tw5(2:end);

log_D_tw5 = log10(D_tw5);
log_rbg_tw5 = log10(rbg_tw5);
P_bgfit_tw5 = polyfit(log_D_tw5,log_rbg_tw5,1);
bgfit_tw5 = P_bgfit_tw5(1)*log_D_tw5 + P_bgfit_tw5(2);
antilog_bgfit_tw5= 10.^(bgfit_tw5);

 
% Bandgap Vs twist plot
%------------

figure1 = figure;
box on

hold on
%scatter(D_tw1,rbg_tw1,'b','Marker','o');
scatter(D_tw2,rbg_tw2,'r','Marker','*');
scatter(D_tw3,rbg_tw3,'g','Marker','sq');
scatter(D_tw4,rbg_tw4,'c','Marker','d');
scatter(D_tw5,rbg_tw5,'m','Marker','+');
%plot(D_tw1,antilog_bgfit_tw1,'b','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw2,antilog_bgfit_tw2,'r','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw3,antilog_bgfit_tw3,'g','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw4,antilog_bgfit_tw4,'c','MarkerSize',10,'LineWidth',1,'LineStyle','-.')
plot(D_tw5,antilog_bgfit_tw5,'m','MarkerSize',10,'LineWidth',1,'LineStyle','-.')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'xscale','log','yscale','log','XTick',[1.25 1.50 1.75 2.00 2.25],'YTick',[0.005 0.01 0.02 0.04 0.08 0.16 0.32],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.3f'))
set(gca,'XTickLabel',{'$1.25$', '$1.50$', '$1.75$', '$2.00$','$2.25$'});
set(gca,'YTickLabel',{'$0.005$', '$0.010$', '$0.020$', '$0.040$', '$0.080$', '$0.160$', '$0.320$'});
xlim([1.24 2.25]);
ylim([0.005 0.32]);
xlabel('$d$ (nm)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\Delta_{Bandgap}$ (eV)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

%ylabel('Relative percentage error (%)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
legend1 = legend('$\alpha = 0.009$','$\alpha = 0.019$','$\alpha = 0.028$','$\alpha = 0.037$','Location','Southeast');
set(legend1,'fontsize',12,'FontName','Times New Roman','Interpreter','LaTex');


set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,'BandgapVsdia','epsc')
hold off
