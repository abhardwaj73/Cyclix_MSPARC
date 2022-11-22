close all
clear
clc

nm2Bohr = 18.89726125;
Ha2eV = 27.211386245988;
indx_ln = 4;

% (13,13) CNT
%----------------

% Radius = 16.6 Bohr, mesh-size=0.15/0.125/0.1, N_c = 13, N_eta = 40, H = 2.309 Bohr

H1 = 2.3091577619; % height of the unit cell
Nc1 = 13;
ndata1 = 6;
rate1 = 0.0004; % rate at which twist is applied in rad/Bohr
tw_d1 = zeros(ndata1,1);
for i = 2:ndata1
    tw_d1(i) = tw_d1(i-1) + rate1 * H1; % twist in rad
end

% Energy at h=0.15 Bohr
% Epl_d1 = [
%           -1.1428424970E+01
%           -1.1428384436E+01
%           -1.1428264167E+01
%           -1.1428063018E+01
%           -1.1427779657E+01 
%           -1.1427413514E+01
%           ] * 2 * Ha2eV * Nc1; % total energy in eV (for ful circle and scaled by 2 to directly obtain torsion modulus from the curve fit)

% Energy at h=0.125 Bohr
% Epl_d1 = [
%           -1.1428410451E+01
%           -1.1428371601E+01
%           -1.1428253238E+01
%           -1.1428053713E+01
%           -1.1427772134E+01
%           -1.1427407939E+01
%          ] * 2 * Ha2eV * Nc1; % 2*total energy in eV

% energy at 0.1 mesh
 Epl_d1 = [
           -1.1428411403E+01
           -1.1428371961E+01
           -1.1428253238E+01
           -1.1428053167E+01
           -1.1427770520E+01
           -1.1427405058E+01
          ] * 2 * Ha2eV * Nc1;

rEpl_d1 = Epl_d1 - Epl_d1(1);

tw_d1 = tw_d1(3:end);
rEpl_d1 = rEpl_d1(3:end);

log_tw_d1 = log10(tw_d1(1:indx_ln));
log_rEpl_d1 = log10(rEpl_d1(1:indx_ln));
P_Efit_d1 = polyfit(log_tw_d1,log_rEpl_d1,1);
Efit_d1 = P_Efit_d1(1)*log_tw_d1 + P_Efit_d1(2);
antilog_Efit_d1= 10.^(Efit_d1);


% (15,0) CNT
%----------------

% Radius = 11 Bohr, mesh-size=0.15/0.125/0.1, N_c = 15, N_eta = 23, H = 4.0017685908 Bohr

H2 = 4.0017685908;
Nc2 = 15;
ndata2 = 6;
rate2 = 0.0004;
tw_d2 = zeros(ndata2,1);
for i = 2:ndata2
    tw_d2(i) = tw_d2(i-1) + rate2 * H2; % input is in rad/Bohr ->rad/nm
end

% energy at 0.15 mesh
% Epl_d2 = [
%           -1.1426121844E+01
%           -1.1426103815E+01 
%           -1.1426048429E+01
%           -1.1425955849E+01
%           -1.1425825900E+01
%           -1.1425658738E+01
% ] * 2 * Ha2eV * Nc2 ; % input is in eV

% energy at 0.125 mesh
% Epl_d2 = [
%           -1.1426083549E+01
%           -1.1426067901E+01
%           -1.1426015106E+01
%           -1.1425925300E+01
%           -1.1425798364E+01 
%           -1.1425634219E+01
% ] * 2 * Ha2eV * Nc2; % input is in eV

%energy at 0.1 mesh
Epl_d2 = [
          -1.1426088912E+01
          -1.1426072903E+01
          -1.1426019688E+01
          -1.1425929262E+01
          -1.1425801588E+01
          -1.1425636659E+01
] * 2 * Ha2eV * Nc2 ; % input is in eV

%4atom unrelaxed initial geometry 
% tw_d2 = tw_d2*2;

% Epl_d2 = [
%            -2.2852041148E+01
%            -2.2852009056E+01
%            -2.2851910901E+01
%            -2.2851742612E+01
%            -2.2851501531E+01
%            -2.2851186264E+01
%  ] * 2 * Ha2eV * Nc2; % input is in eV

% 2-atom unrelaxed initial geometry
% Epl_d2 = [
%            -1.1426069295E+01
%            -1.1426054674E+01
%            -1.1426002918E+01
%            -1.1425914110E+01
%            -1.1425788066E+01
%            -1.1425624700E+01
%   ] * 2 * Ha2eV * Nc2; % input is in eV


rEpl_d2 = Epl_d2 - Epl_d2(1);

tw_d2 = tw_d2(3:end);
rEpl_d2 = rEpl_d2(3:end);

log_tw_d2 = log10(tw_d2(1:indx_ln));
log_rEpl_d2 = log10(rEpl_d2(1:indx_ln));
P_Efit_d2 = polyfit(log_tw_d2,log_rEpl_d2,1);
Efit_d2 = P_Efit_d2(1)*log_tw_d2 + P_Efit_d2(2);
antilog_Efit_d2= 10.^(Efit_d2);

% (20,0) CNT
%----------------

% Radius = 14.7 Bohr, mesh-size=0.15/0.125/0.1, N_c = 20, N_eta = 23, H = 4.0002688791 Bohr

H3 = 4.0002688791;
Nc3 = 20;
ndata3 = 6;
rate3 = 0.0004;
tw_d3 = zeros(ndata3,1);
for i = 2:ndata3
    tw_d3(i) = tw_d3(i-1) + rate3 * H3; % input is in rad/Bohr ->rad/nm
end

% energy at 0.15 mesh
% Epl_d3 = [
%           -1.1427970567E+01
%           -1.1427934417E+01
%           -1.1427833151E+01
%           -1.1427667038E+01
%           -1.1427436071E+01
%           -1.1427140771E+01
% ] * 2 * Ha2eV * Nc3; % input is in eV

% energy at 0.125 mesh
% Epl_d3 = [
%            -1.1427933283E+01
%            -1.1427899998E+01
%            -1.1427802124E+01
%            -1.1427639808E+01
%            -1.1427412885E+01
%            -1.1427121295E+01
% ] * 2 * Ha2eV * Nc3; % input is in eV

% energy at 0.1 mesh
Epl_d3 = [
          -1.1427937009E+01 
          -1.1427903884E+01
          -1.1427805898E+01
          -1.1427643004E+01
          -1.1427415249E+01
          -1.1427122704E+01
] * 2 * Ha2eV * Nc3; % input is in eV

%4atom unrelaxed initial geometry
% tw_d3 = tw_d3*2;

% Epl_d3 = [
%            -2.2855828574E+01 
%            -2.2855764698E+01 
%            -2.2855572654E+01 
%            -2.2855251999E+01 
%            -2.2854801932E+01
%            -2.2854221589E+01
%  ] * 2 * Ha2eV * Nc3; % input is in eV

% 2atom initial unrelaxed geometry
% Epl_d3 = [
%            -1.1427928561E+01 
%            -1.1427895210E+01
%            -1.1427797104E+01
%            -1.1427634398E+01
%            -1.1427406888E+01
%            -1.1427114583E+01 
%  ] * 2 * Ha2eV * Nc3; % input is in eV

rEpl_d3 = Epl_d3 - Epl_d3(1);

tw_d3 = tw_d3(3:end);
rEpl_d3 = rEpl_d3(3:end);

log_tw_d3 = log10(tw_d3(1:indx_ln));
log_rEpl_d3 = log10(rEpl_d3(1:indx_ln));
P_Efit_d3 = polyfit(log_tw_d3,log_rEpl_d3,1);
Efit_d3 = P_Efit_d3(1)*log_tw_d3 + P_Efit_d3(2);
antilog_Efit_d3= 10.^(Efit_d3);

% (16,0) CNT
%----------------

% Radius = 11.77 Bohr, mesh-size=0.1, N_c = 16, N_eta = 23, H = 4.0023317752 Bohr

H6 = 4.0023317752;
Nc6 = 16;
ndata6 = 6;
rate6 = 0.0004;
tw_d6 = zeros(ndata6,1);
for i = 2:ndata6
    tw_d6(i) = tw_d6(i-1) + rate6 * H6; % input is in rad/Bohr ->rad/nm
end

% energy at 0.1 mesh
Epl_d6 = [
          -1.1426624938E+01
          -1.1426605087E+01
          -1.1426544108E+01
          -1.1426441816E+01
          -1.1426298131E+01
          -1.1426112855E+01
          %-1.1425885763E+01
] * 2 * Ha2eV * Nc6; % input is in eV

rEpl_d6 = Epl_d6 - Epl_d6(1);

tw_d6 = tw_d6(3:end);
rEpl_d6 = rEpl_d6(3:end);

log_tw_d6 = log10(tw_d6(1:indx_ln));
log_rEpl_d6 = log10(rEpl_d6(1:indx_ln));
P_Efit_d6 = polyfit(log_tw_d6,log_rEpl_d6,1);
Efit_d6 = P_Efit_d6(1)*log_tw_d6 + P_Efit_d6(2);
antilog_Efit_d6= 10.^(Efit_d6);

% (27,0) CNT
%----------------

% Radius = 19.84 Bohr, mesh-size=0.15, N_c = 27, N_eta = 23, H = 4.0019793139 Bohr

H7 = 4.0019793139;
Nc7 = 27;
ndata7 = 6;
rate7 = 0.0004;
tw_d7 = zeros(ndata7,1);
for i = 2:ndata7
    tw_d7(i) = tw_d7(i-1) + rate7 * H7; % input is in rad/Bohr ->rad/nm
end

% energy at 0.1 mesh
Epl_d7 = [
          -1.1428997073E+01
          -1.1428939510E+01
          -1.1428761891E+01
          -1.1428464173E+01
          -1.1428046431E+01
          -1.1427508876E+01
          %-1.1426851724E+01 
] * 2 * Ha2eV * Nc7; % input is in eV

rEpl_d7 = Epl_d7 - Epl_d7(1);

tw_d7 = tw_d7(3:end);
rEpl_d7 = rEpl_d7(3:end);

log_tw_d7 = log10(tw_d7(1:indx_ln));
log_rEpl_d7 = log10(rEpl_d7(1:indx_ln));
P_Efit_d7 = polyfit(log_tw_d7,log_rEpl_d7,1);
Efit_d7 = P_Efit_d7(1)*log_tw_d7 + P_Efit_d7(2);
antilog_Efit_d7= 10.^(Efit_d7);

% Twist energy Vs twist plot
%------------

figure1 = figure;
box on

hold on
scatter(tw_d2,rEpl_d2,'k','Marker','o');
scatter(tw_d6,rEpl_d6,'k','Marker','+');
scatter(tw_d5,rEpl_d5,'k','Marker','^');
scatter(tw_d3,rEpl_d3,'k','Marker','*');
scatter(tw_d4,rEpl_d4,'k','Marker','x');
scatter(tw_d1,rEpl_d1,'k','Marker','d');
scatter(tw_d7,rEpl_d7,'k','Marker','.');


plot(tw_d2,antilog_Efit_d2,'k','Markersize',10,'LineWidth',1,'LineStyle','-');
plot(tw_d6,antilog_Efit_d6,'k','Markersize',10,'LineWidth',1,'LineStyle','-');
plot(tw_d5,antilog_Efit_d5,'k','Markersize',10,'LineWidth',1,'LineStyle','-');
plot(tw_d3,antilog_Efit_d3,'k','Markersize',10,'LineWidth',1,'LineStyle','-');
plot(tw_d4,antilog_Efit_d4,'k','Markersize',10,'LineWidth',1,'LineStyle','-');
plot(tw_d1,antilog_Efit_d1,'k','Markersize',10,'LineWidth',1,'LineStyle','-');
plot(tw_d7,antilog_Efit_d7,'k','Markersize',10,'LineWidth',1,'LineStyle','-');

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'xscale','log','yscale','log','XTick',[0.001 0.002 0.004 0.008],'XMinorTick','on','YTick',[1e-2 1e-1 1e0],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.3f'))
%set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
set(gca,'XTickLabel',{'$0.001$', '$0.002$', '$0.004$','$0.008$'});
set(gca,'YTickLabel',{'$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
xlim([0.0009 0.01]);
ylim([1e-2 1e0]);

xlabel('$\alpha$ (rad)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$2\mathcal{E}_{twist}$ (eV)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

%ylabel('Relative percentage error (%)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
legend1 = legend('$d = 1.16$ nm','$d = 1.24$ nm','$d = 1.41$ nm','$d = 1.55$ nm','$d = 1.61$ nm','$d = 1.75$ nm','$d = 2.10$ nm','Location','Northwest');
set(legend1,'fontsize',14,'FontName','Times New Roman','Interpreter','LaTex');

%str = 'Slope = 2.0';
%dim = [.2 .58 .3 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'twistenergyVstwist_linear','epsc')
hold off
