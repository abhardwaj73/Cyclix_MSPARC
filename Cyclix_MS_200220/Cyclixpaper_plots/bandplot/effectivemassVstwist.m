close all
clear
clc

nm2Bohr = 18.89726125;
Ha2eV = 27.211386245988;

%Chiral angle
omega = @(n,m) acos((2*n+m)/(2*sqrt(n^2+m^2+n*m)));

%===============================================
% Tw_n = Tw_(n-1) + (n-1) * rate ; n = 1,2, ...
%===============================================

% (13,13) CNT
%----------------

H1 = [2.3097769325]/nm2Bohr;
rad1 = 0.875;

omega1 = omega(13,13);

tw_d1 = [2.4166097335E-01; 2.4258525219E-01;2.4350893866E-01;2.4443309635E-01;2.4535720473E-01;2.4628142786E-01]/H1(1);

tw_d1 = tw_d1 - tw_d1(1);

shr_st1 = tw_d1*rad1*100;

me_d1 = [0.0034;0.0082;0.0173;0.0263;0.0346;0.0442];

mh_d1 = -[0.0034;0.0082;0.0172;0.0259;0.0342;0.0438];

ve_d1 = [8.66*1e5];

vh_d1 = [8.95*1e5];

bg_d1 = [0.0006497552;0.0025793338;0.0050977869;0.0075473631;0.0099754528;0.0125705660] * Ha2eV;

[Pm1,Sm1] = polyfit(shr_st1,(me_d1-me_d1(1)),1);
[Pg1,Sg1] = polyfit(shr_st1,(bg_d1-bg_d1(1)),1);

% (20,0) CNT
%----------------

H2 = [4.0012612687]/nm2Bohr;
rad2 = 0.78;

omega2 = omega(20,0);

tw_d2 = [1.5707963268E-01;1.5867974107E-01;1.6027993216E-01;1.6188023681E-01;1.6348050953E-01;1.6508109060E-01]/H2(1);

tw_d2 = tw_d2 - tw_d2(1);

shr_st2 = tw_d2*rad2*100;

me_d2 = [0.0743;0.0772;0.0761;0.0709;0.0741;0.0730];

mh_d2 = -[0.0719;0.0745;0.0734;0.0685;0.0717;0.0711];

ve_d2 = [0];

vh_d2 = [0];

bg_d2 = [0.0177922114;0.0183717523;0.0181456737;0.0171327177;0.0177679752;0.0176048954] * Ha2eV;

Pm2 = polyfit(shr_st2,(me_d2-me_d2(1)),1);
Pg2 = polyfit(shr_st2,(bg_d2-bg_d2(1)),1);

% (16,0) CNT
%----------------

H3 = [4.0023317752]/nm2Bohr;
rad3 = 0.625;

omega3 = omega(16,0);

tw_d3 = [1.9634954085E-01;1.9795047356E-01;1.9955133873E-01;2.0115242157E-01;2.0275347201E-01;2.0435463346E-01]/H3(1);

tw_d3 = tw_d3 - tw_d3(1);

shr_st3 = tw_d3*rad3*100;

me_d3 = [0.0593;0.0590;0.0595;0.0594;0.0598;0.0599];

mh_d3 = -[0.0586;0.0583;0.0588;0.0585;0.0590;0.0593];

ve_d3 = [0];

vh_d3 = [0];

bg_d3 = [0.0205869485;0.0204728694;0.0207437207;0.0206586890;0.0207721590;0.0208740889] * Ha2eV;

Pm3 = polyfit(shr_st3,(me_d3-me_d3(1)),1);
Pg3 = polyfit(shr_st3,(bg_d3-bg_d3(1)),1);

% (27,0) CNT
%----------------

H4 = [4.0019793139]/nm2Bohr;
rad4 = 1.05;

omega4 = omega(27,0);

tw_d4 = [1.1635528347E-01;1.1795607519E-01;1.1955773183E-01;1.2115896609E-01;1.2276071629E-01]/H4(1);%;1.2436320551E-01

tw_d4 = tw_d4 - tw_d4(1);

shr_st4 = tw_d4*rad4*100;

me_d4 = [0.0041;0.0052;0.0050;0.0066;0.0077];%;0.0077

mh_d4 = -[0.0041;0.0052;0.0050;0.0065;0.0077];%;0.0077

ve_d4 = [8.34*1e5];

vh_d4 = [8.47*1e5];

bg_d4 = [0.0009087744;0.0016503023;0.0012836314;0.0017223591;0.0025018397] * Ha2eV;%;0.0019300801

Pm4 = polyfit(shr_st4,(me_d4-me_d4(1)),1);
Pg4 = polyfit(shr_st4,(bg_d4-bg_d4(1)),1);

% (22,11) CNT
%-------------

H5 = [1.5124912124]/nm2Bohr;
rad5 = 1.13;

omega5 = omega(22,11);

tw_d5 = [2.0399952296E-01;2.0460451945E-01;2.04907117E-01;2.0520955054E-01]/H5(1);%;2.0581497183E-01;2.0642016346E-01;2.0702566820E-01

tw_d5 = tw_d5 - tw_d5(1);

shr_st5 = tw_d5*rad5*100;

me_d5 = [0.0433;0.0540;0.0606;0.0622];%;0.0510;0.0450];%;0.0369

mh_d5 = -[0.0432;0.0530;0.0594;0.0612];%;0.0504;0.0450];%;0.0364

ve_d5 = [0];

vh_d5 = [0];

bg_d5 = [0.0118583766;0.0144859116;0.0160;0.0164725983]* Ha2eV;%0.0144009442;0.0115158506;0.0090462582] * Ha2eV; 

Pm5 = polyfit(shr_st5,(me_d5-me_d5(1)),1);
Pg5 = polyfit(shr_st5,(bg_d5-bg_d5(1)),1);

% (18,9) CNT
%-----------

H6 = [1.5130233554]/nm2Bohr;
rad6 = 0.925;

omega6 = omega(18,9);

tw_d6 = [2.4933275028E-01;2.4993795963E-01;2.5054291064E-01;2.5114786009E-01;2.5175292042E-01]/H6(1);%;2.5235795549E-01

tw_d6 = tw_d6 - tw_d6(1);

shr_st6 = tw_d6*rad6*100;

me_d6 = [0.0051;0.0065;0.0150;0.0230;0.0324];%;0.0389

mh_d6 = -[0.0051;0.0065;0.0150;0.0228;0.0322];%;0.0386

ve_d6 = [8.16*1e5];

vh_d6 = [8.10*1e5];

bg_d6 = [0.0011349332;0.0020306156;0.0043571108;0.0063898496;0.0089721109] * Ha2eV;%;

Pm6 = polyfit(shr_st6,(me_d6-me_d6(1)),1);
Pg6 = polyfit(shr_st6,(bg_d6-bg_d6(1)),1);

% Effective mass Vs twist plot
%-----------------------------

figure1 = figure;
box on

hold on
% S1 = plot(shr_st3,me_d3,'d','Marker','o','Markersize',12);S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560];S1.MarkerFaceColor = [0.4940, 0.1840, 0.5560];
% S2 = plot(shr_st2,me_d2,'d','Marker','p','Markersize',12);S2.MarkerEdgeColor = [0 0.6 0];S2.MarkerFaceColor = [0 0.6 0];
% S3 = plot(shr_st1(2:end),me_d1(2:end),'d','Marker','d','Markersize',12);S3.MarkerEdgeColor = [0 0 1];S3.MarkerFaceColor = [0 0 1];
% S4 = plot(shr_st6(2:end),me_d6(2:end),'d','Marker','h','Markersize',12);S4.MarkerEdgeColor = [1 0 0];S4.MarkerFaceColor = [1 0 0];
% S5 = plot(shr_st4(2:end),me_d4(2:end),'d','Marker','sq','Markersize',12);S5.MarkerEdgeColor = [0.64 0.08 0.18];S5.MarkerFaceColor = [0.64 0.08 0.18];
% S6 = plot(shr_st5,me_d5,'d','Marker','^','Markersize',12);S6.MarkerEdgeColor = [0.95 0.00 0.75];S6.MarkerFaceColor = [0.95 0.00 0.75];

S1 = plot(shr_st3,abs(mh_d3),'d','Marker','o','Markersize',12);S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560];
S2 = plot(shr_st2,abs(me_d2),'d','Marker','p','Markersize',12);S2.MarkerEdgeColor = [0 0.6 0];
S3 = plot(shr_st1(2:end),abs(mh_d1(2:end)),'d','Marker','d','Markersize',12);S3.MarkerEdgeColor = [0 0 1];
S4 = plot(shr_st6(2:end),abs(mh_d6(2:end)),'d','Marker','h','Markersize',12);S4.MarkerEdgeColor = [1 0 0];
S5 = plot(shr_st4(2:end),abs(mh_d4(2:end)),'d','Marker','sq','Markersize',12);S5.MarkerEdgeColor = [0.64 0.08 0.18];
S6 = plot(shr_st5,abs(mh_d5),'d','Marker','^','Markersize',12);S6.MarkerEdgeColor = [0.95 0.00 0.75];

% Law line

c1 = (Pm1(1)/sin(3*omega1)+Pm5(1)/sin(3*omega5)+Pm6(1)/sin(3*omega6))/3*100
m = @ (m0,st,a) m0 + c1*(st/100)*sin(3*a);

plot(shr_st3,m(me_d3(1),shr_st3,omega3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',1,'LineStyle','-')
plot(shr_st2,m(me_d2(1),shr_st2,omega2),'color',[0 0.6 0],'LineWidth',1,'LineStyle','-')
plot(shr_st1,m(0,shr_st1,omega1),'color',[0 0 1],'LineWidth',1,'LineStyle','-')
plot(shr_st6,m(0,shr_st6,omega6),'color',[1 0 0],'LineWidth',1,'LineStyle','-')
plot(shr_st4,m(me_d4(2),shr_st4,omega4),'color',[0.64 0.08 0.18],'LineWidth',1,'LineStyle','-')
plot(shr_st5,m(me_d5(1),shr_st5,omega5),'color',[0.95 0.00 0.75],'LineWidth',1,'LineStyle','-')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'XTick',[0 1 2 3],'XMinorTick','on','YTick',[0.00 0.02 0.04 0.06 0.08],'YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',24);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.2f'))
%set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
%set(gca,'YTickLabel',{'$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
xlim([0 3.4]);
ylim([0.00 0.09]);

xlabel('Shear strain (\%)','FontSize',24,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Effective mass of charge carriers ($m_e$)','FontSize',24,'FontName','Times New Roman','Interpreter','LaTex');



% yyaxis right

% plot(tw_d3,bg_d3,'color',[0.4940, 0.1840, 0.5560],'Marker','o','MarkerFaceColor',[0.4940, 0.1840, 0.5560],'Markersize',12,'LineWidth',1,'LineStyle','-.')
% plot(tw_d2,bg_d2,'color',[0 0.6 0],'Marker','p','MarkerFaceColor',[0 0.6 0],'Markersize',12,'LineWidth',1,'LineStyle','-.')
% plot(tw_d1,bg_d1,'color',[0 0 1],'Marker','d','MarkerFaceColor',[0 0 1],'Markersize',12,'LineWidth',1,'LineStyle','-.')
% plot(tw_d6,bg_d6,'color',[1 0 0],'Marker','h','MarkerFaceColor',[1 0 0],'Markersize',12,'LineWidth',1,'LineStyle','-.')
% plot(tw_d4,bg_d4,'color',[0.64 0.08 0.18],'Marker','sq','MarkerFaceColor',[0.64, 0.08, 0.18],'Markersize',12,'LineWidth',1,'LineStyle','-.')
% plot(tw_d5,bg_d5,'color',[0.95 0.00 0.75],'Marker','^','MarkerFaceColor',[0.95 0.00 0.75],'Markersize',12,'LineWidth',1,'LineStyle','-.')

% set(gca,'TickLabelInterpreter','LaTex');
% set(gca,'Ycolor',[0 0 0],'YTick',[0.0 0.1 0.2 0.3 0.4 0.5 0.6],'YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
% set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
% ylim([0.0 0.65]);

% ylabel('Bandgap (eV)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');






legend1 = legend('$(16,0)$','$(20,0)$','$(13,13)$','$(18,9)$','$(27,0)$','$(22,11)$','Location','North','Orientation','Horizontal');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

%str = 'Slope = 2.0';
%dim = [.2 .58 .3 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPosition',[pos(1) pos(2) pos(3)*1.4 pos(4)*1.2],'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'effectivemassVstwist_small','epsc')
hold off


figure2 = figure;
box on

hold on

S1 = plot(shr_st3,bg_d3,'d','Marker','o','Markersize',12);S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560];
S2 = plot(shr_st2,bg_d2,'d','Marker','p','Markersize',12);S2.MarkerEdgeColor = [0 0.6 0];
S3 = plot(shr_st1,bg_d1,'d','Marker','d','Markersize',12);S3.MarkerEdgeColor = [0 0 1];
S4 = plot(shr_st6,bg_d6,'d','Marker','h','Markersize',12);S4.MarkerEdgeColor = [1 0 0];
S5 = plot(shr_st4,bg_d4,'d','Marker','sq','Markersize',12);S5.MarkerEdgeColor = [0.64 0.08 0.18];
S6 = plot(shr_st5,bg_d5,'d','Marker','^','Markersize',12);S6.MarkerEdgeColor = [0.95 0.00 0.75];

t0 = (Pg1(1)/sin(3*omega1)+Pg5(1)/sin(3*omega5)+Pg6(1)/sin(3*omega6))/3/3*100
bg = @ (bg0,st,a) bg0 + 3*t0*(st/100)*sin(3*a);

plot(shr_st3,bg(bg_d3(1),shr_st3,omega3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',1,'LineStyle','-')
plot(shr_st2,bg(bg_d2(1),shr_st2,omega2),'color',[0 0.6 0],'LineWidth',1,'LineStyle','-')
plot(shr_st1,bg(bg_d1(1),shr_st1,omega1),'color',[0 0 1],'LineWidth',1,'LineStyle','-')
plot(shr_st6,bg(0.01,shr_st6,omega6),'color',[1 0 0],'LineWidth',1,'LineStyle','-')
plot(shr_st4,bg(bg_d4(1),shr_st4,omega4),'color',[0.64 0.08 0.18],'LineWidth',1,'LineStyle','-')
plot(shr_st5,bg(bg_d5(1),shr_st5,omega5),'color',[0.95 0.00 0.75],'LineWidth',1,'LineStyle','-')


set(gca,'TickLabelInterpreter','LaTex');
set(gca,'Ycolor',[0 0 0],'YTick',[0.0 0.1 0.2 0.3 0.4 0.5 0.6],'YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',24);
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
ylim([0.0 0.65]);



set(gca,'XTick',[0 1 2 3],'XMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',24);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
xlim([0 3.4]);

xlabel('Shear strain (\%)','FontSize',24,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Bandgap (eV)','FontSize',24,'FontName','Times New Roman','Interpreter','LaTex');


legend2 = legend('$(16,0)$','$(20,0)$','$(13,13)$','$(18,9)$','$(27,0)$','$(22,11)$','Location','North','Orientation','Horizontal');
set(legend2,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

%str = 'Slope = 2.0';
%dim = [.2 .58 .3 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure2,'Units','Inches');
pos = get(figure2,'Position');
set(figure2,'PaperPosition',[pos(1) pos(2) pos(3)*1.4 pos(4)*1.2],'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'bandgapVstwist_small','epsc')
hold off
