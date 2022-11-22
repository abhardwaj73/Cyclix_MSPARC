close all
clear
clc

nm2Bohr = 18.89726125;
Ha2eV = 27.211386245988;
%indx_ln = 5;

% (18,9) CNT
%----------------

% Radius = 17.46 Bohr, mesh-size=0.15, N_c = 9, N_eta = 70, H = 1.5 Bohr

H = [1.5130233554;1.5123736013;1.5126508325;1.5125876811;1.5130391119;1.5132178566;1.5135765153]/nm2Bohr;%1.5123775448;
Nc = 9;
rate = 0.0004;
tw = [
         2.4933275028E-01
         %2.4993795963E-01
         2.5054291064E-01
         2.5114786009E-01
         2.5175292042E-01
         2.5235795549E-01
         ]/H(1);

% energy at 0.15 mesh
E1 = [
          -1.1428547731E+01
          %-1.1428555894E+01
          -1.1428431724E+01
          -1.1428210544E+01
          -1.1427925651E+01
          -1.1427527839E+01
] * 2 * Ha2eV * Nc/H(1) ; % input is in Ha -> eV

bg1 = [0.0011349332;0.0020306156;0.0043571108;0.0063898496;0.0089721109] * Ha2eV;

 % without stress relaxation
 E2 = [
        -1.1428613564E+01
        %-1.1428568211E+01
        -1.1428438702E+01
        -1.1428224436E+01
        -1.1427924479E+01
        -1.1427536678E+01
 ] * 2 * Ha2eV * Nc/H(1) ;

bg2 = [0.0053;0.0044;0.00496;0.006691;0.00868;0.01093]
rE1 = E1 - E1(1);
tw = tw - tw(1);
tw = tw(2:end);
tw2 = tw.^2;
rE1 = rE1(2:end);

log_tw = log10(tw);
log_rE1 = log10(rE1);
P_Efit1 = polyfit(log_tw,log_rE1,1);
Efit1 = P_Efit1(1)*log_tw + P_Efit1(2);
antilog_Efit1= 10.^(Efit1);


rE2 = E2 - E2(1);
rE2 = rE2(2:end);

log_rE2 = log10(rE2);
P_Efit2 = polyfit(log_tw,log_rE2,1);
Efit2 = P_Efit2(1)*log_tw + P_Efit2(2);
antilog_Efit2= 10.^(Efit2);

% Twist energy Vs twist plot
%------------

figure1 = figure;
box on

hold on
%scatter(tw_d2,rEpl_d2,'Marker','o');
plot(tw,rE1*0.5,'r','Marker','o','MarkerSize',12);
plot(tw,rE2*0.5,'b','Marker','*','MarkerSize',12); 

%plot(tw_d2,antilog_Efit_d2,'LineWidth',1,'LineStyle','-');
%plot(tw,antilog_Efit1*0.5,'color','r','LineWidth',1,'LineStyle','-');
%plot(tw,antilog_Efit2*0.5,'color','b','LineWidth',1,'LineStyle','-');

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'XTick',[0.00 0.01 0.02 0.03 0.04],'YTick',[0 1.0 2.0 3.0 4.0],'XMinorTick','on','YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'));
%set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
set(gca,'XTickLabel',{'$0.00$', '$0.01$', '$0.02$','$0.03$','$0.04$'});
xlim([0.00 0.04]);
ylim([0 4]);

xlabel('Twist (rad/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\mathcal{E}_{twist}$ (eV/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

legend1 = legend('with stress relax.','w/o stress relax.','Location','Northwest');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

%str = 'Slope = 2.0';
%dim = [.2 .58 .3 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4 2]);
saveas(gcf,'effect_stress_relax','epsc');
hold off



