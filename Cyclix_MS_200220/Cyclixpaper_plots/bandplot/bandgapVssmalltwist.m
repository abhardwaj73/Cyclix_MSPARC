close all
clear
clc

nm2Bohr = 18.89726125;
Ha2eV = 27.211386245988;
%===============================================
% Tw_n = Tw_(n-1) + (n-1) * rate ; n = 1,2, ...
%===============================================

% (13,13) CNT
%----------------

% Radius = 16.6 Bohr, mesh-size=0.3, N_c = 13, N_eta = 40, H = 2.309 Bohr
H1 = [2.3097769325]/nm2Bohr;%;2.4326444319

Nc1 = 13;
ndata1 = 6;

tw_d1 = [2.4166097335E-01; 2.4258525219E-01;2.4350893866E-01;2.4443309635E-01;2.4535720473E-01;2.4628142786E-01]/H1(1);%;2.9328548535E-01

tw_d1 = tw_d1 - tw_d1(1);

bg_d1 = [0.0006497552;0.0025793338;0.0050977869;0.0075473631;0.0099754528;0.0125705660] * Ha2eV;

kpt_homo = [0.000000 -0.3350];

kpt_lumo = [0.000000 -0.3350];

% (20,0) CNT
%----------------

% Radius = 14.7 Bohr, mesh-size=0.3, N_c = 20, N_eta = 23, H = 4.0002688791 Bohr

H2 = [4.0012612687]/nm2Bohr;

Nc2 = 20;
ndata2 = 7;

tw_d2 = [1.5707963268E-01;1.5867974107E-01;1.6027993216E-01;1.6188023681E-01;1.6348050953E-01;1.6508109060E-01]/H2(1);

tw_d2 = tw_d2 - tw_d2(1);

bg_d2 = [0.0177922114;0.0183717523;0.0181456737;0.0171327177;0.0177679752;0.0176048954] * Ha2eV;

kpt_homo = [0.350000 -0.325652];
kpt_lumo = [0.350000 -0.325652];
% (16,0) CNT
%----------------

% Radius = 11.77 Bohr, mesh-size=0.3, N_c = 16, N_eta = 23, H = 4.0007020478 Bohr

H3 = [4.0023317752]/nm2Bohr;

Nc3 = 16;
ndata3 = 7;
tw_d3 = [1.9634954085E-01;1.9795047356E-01;1.9955133873E-01;2.0115242157E-01;2.0275347201E-01;2.0435463346E-01]/H3(1);

tw_d3 = tw_d3 - tw_d3(1);

% bg_d3 = [0.0243857345;0.0222633268;0.0212922933;0.0260218015;0.0338100413;
%          0.0223978280;0.0143157350] * Ha2eV; % input is in eV

% kpt_homo = [0.312500 -0.333333;0.312500 -0.333333;0.312500 -0.333333;0.312500 -0.333333;0.625000 0.266667;
%             0.375000 -0.233333;0.375000 -0.233333];

% kpt_lumo = [0.312500 -0.333333;0.312500 -0.333333;0.312500 -0.333333;0.312500 -0.333333;0.625000 0.266667;
%             0.375000 -0.233333;0.375000 -0.233333];

bg_d3 = [0.0205869485;0.0204728694;0.0207437207;0.0206586890;0.0207721590;0.0208740889] * Ha2eV; % input is in eV

kpt_homo = [0.312500 -0.343043];

kpt_lumo = [0.312500 -0.343043];

% (27,0) CNT
%----------------

% Radius = 19.84 Bohr, mesh-size=0.3, N_c = 27, N_eta = 23, H = 4.0025688504 Bohr

H4 = [4.0019793139]/nm2Bohr;

Nc4 = 27;
ndata4 = 6;

tw_d4 = [1.1635528347E-01;1.1795607519E-01;1.1955773183E-01;1.2115896609E-01;1.2276071629E-01;1.2436320551E-01]/H4(1);

tw_d4 = tw_d4 - tw_d4(1);

bg_d4 = [0.0009087744;0.0016503023;0.0012836314;0.0017223591;0.0025018397;0.0019300801] * Ha2eV;

kpt_homo = [0.333333 -0.332609];

kpt_lumo = [0.333333 -0.332609];

% (22,11) CNT
%-------------

% Radius = 21.4 Bohr, mesh-size=0.3, N_c = 11, N_eta = 30, H = 1.51 Bohr

H5 = [1.5124912124]/nm2Bohr;

Nc5 = 11;
ndata5 = 7;

tw_d5 = [2.0399952296E-01;2.0460451945E-01;2.0520955054E-01;2.0581497183E-01;2.0642016346E-01;2.0702566820E-01]/H5(1);

tw_d5 = tw_d5 - tw_d5(1);

%bg_d5 = [0.0393185595;0.0434604405;0.0503994408;0.0464914867;0.0318319969;
%         0.0221928810;0.0073213476] * Ha2eV;

bg_d5 = [0.0118583766;0.0144859116;0.0164725983;0.0144009442;0.0115158506;0.0090462582] * Ha2eV;         

kpt_homo = [0.363636 0.344857];

kpt_lumo = [0.363636 0.344857];

% (18,9) CNT
%-----------

% Radius = 17.46 Bohr, mesh-size=0.3, N_c = 9, N_eta = 30, H = 1.51 Bohr

H6 = [1.5130233554]/nm2Bohr;%1.5547457030;

Nc6 = 9;
ndata6 = 7;

tw_d6 = [2.4933275028E-01;2.4993795963E-01;2.5054291064E-01;2.5114786009E-01;2.5175292042E-01]/H6(1);

tw_d6 = tw_d6 - tw_d6(1);

% bg_d6 = [0.0053197103;0.0069498411;0.0133157873;0.0245182875;0.0163234359;
%          0.0175134516;0.0306393328;0.0254754378] * Ha2eV;

bg_d6 = [0.0011349332;0.0020306156;0.0043571108;0.0063898496;0.0089721109] * Ha2eV;%0.0115604264;         

kpt_homo = [0.333333 0.334571];

kpt_lumo = [0.333333 0.334571];

% Band gap Vs twist plot
%------------

figure1 = figure;
box on

hold on
plot(tw_d3,bg_d3,'color',[0.4940, 0.1840, 0.5560],'Marker','o','Markersize',12,'LineWidth',1,'LineStyle','-')
plot(tw_d2,bg_d2,'color',[0 0.6 0],'Marker','p','Markersize',12,'LineWidth',1,'LineStyle','-')
plot(tw_d1,bg_d1,'color',[0 0 1],'Marker','d','Markersize',12,'LineWidth',1,'LineStyle','-')
plot(tw_d6,bg_d6,'color',[1 0 0],'Marker','h','Markersize',12,'LineWidth',1,'LineStyle','-')
plot(tw_d4,bg_d4,'color',[0.64 0.08 0.18],'Marker','sq','Markersize',12,'LineWidth',1,'LineStyle','-')
plot(tw_d5,bg_d5,'color',[0.95 0.00 0.75],'Marker','^','Markersize',12,'LineWidth',1,'LineStyle','-')


set(gca,'TickLabelInterpreter','LaTex');
set(gca,'XTick',[0.00 0.01 0.02 0.03 0.04],'XMinorTick','on','YTick',[0.0 0.2 0.4 0.6],'YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
%set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
%set(gca,'YTickLabel',{'$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
xlim([0.00 0.04]);
ylim([0.0 0.65]);

xlabel('Twist (rad/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Bandgap (eV)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

legend1 = legend('$(16,0)$','$(20,0)$','$(13,13)$','$(18,9)$','$(27,0)$','$(22,11)$','Location','North','Orientation','horizontal');
set(legend1,'fontsize',11,'FontName','Times New Roman','Interpreter','LaTex');

%str = 'Slope = 2.0';
%dim = [.2 .58 .3 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'bandgapVssmalltwist','epsc')
hold off
