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

% Radius = 16.6 Bohr, mesh-size=0.15, N_c = 13, N_eta = 40, H = 2.309 Bohr

H1 = [2.309;2.310;2.310;2.3101092988;2.3105697652;
      2.3117082998;2.3143331865;2.3152896439;2.3175341149;2.3223812442;
      2.3264483065;2.3285365458;2.3335400767;2.3424470641;2.3565069034;
      2.3610363182;2.3699167072;2.3871009369;2.4358926927;2.4466551027;
      2.4609627568;3.9901396451]/nm2Bohr;

st1 = (H1 - H1(1))*100/H1(1);

R1 = [16.597;16.594;16.598;16.624;16.624;16.625;16.620;16.620;16.620;16.620;16.617;16.617;16.617;16.618;16.654;16.632;16.689;16.738;16.853;16.886;16.701];
Nc1 = 13;
ndata1 = 22;
rate1 = 0.0004;

tw_d1 = [2.4166097335E-01;2.4397291607E-01;2.4512871480E-01;2.4859637249E-01;2.4975196922E-01;
         2.5206317343E-01;2.5553097575E-01;2.5668767554E-01;2.5900174059E-01;2.6247540102E-01;
         2.6478619331E-01;2.6594713960E-01;2.6827077114E-01;2.7176170465E-01;2.7639636730E-01;
         2.7756887363E-01;2.7991735886E-01;2.8345199407E-01;2.9046671770E-01;2.9167229519E-01;
         2.9280000000E-01;2.9400000000E-01]/H1(1);

tw_d1 = tw_d1 - tw_d1(1);

Epl_d1 = [-1.1428354191E+01;-1.1428117326E+01;-1.1427821601E+01;-1.1426228953E+01;-1.1425447593E+01;
          -1.1423519790E+01;-1.1419721061E+01;-1.1418212407E+01;-1.1414839603E+01;-1.1408907377E+01;
          -1.1404388558E+01;-1.1401942403E+01;-1.1396696702E+01;-1.1387962324E+01;-1.1374925236E+01;
          -1.1371319955E+01;-1.1363967525E+01;-1.1352067302E+01;-1.1326772603E+01;-1.1322162870E+01;
          -1.1318084424E+01;-1.1394123652E+01] * Ha2eV * Nc1/H1(1);

rEpl_d1 = Epl_d1 - Epl_d1(1);


ppval(fnder(spline(tw_d1,rEpl_d1),2),tw_d1);%17



% (20,0) CNT
%----------------

% Radius = 14.7 Bohr, mesh-size=0.15, N_c = 20, N_eta = 23, H = 4.0002688791 Bohr

H3 = [3.999;4.001;4.0032281599;4.004;4.0041619257;
      4.0054856520;4.0084784881;4.0165531515;4.0188789056;4.0247863602;
      4.0354123302;4.0436899184;4.0477129326;4.0570970967;4.0736633130;
      4.0977345770;4.1044896317;4.1175675337;4.1393605870;4.1912550617;
      4.2042904672;4.3394605718;4.7992594608;]/nm2Bohr;

st3 = (H3 - H3(1))*100/H3(1);

Nc3 = 20;
ndata3 = 11;
rate3 = 0.0004;

tw_d3 = [1.5707963268E-01;1.6108084942E-01;1.6308237428E-01;1.6708585683E-01;1.6908911484E-01;
         1.7109155663E-01;1.7509777500E-01;1.8110477253E-01;1.8311329764E-01;1.8713272246E-01;
         1.9316994823E-01;1.9717319825E-01;1.9919459538E-01;2.0324173615E-01;2.0932763406E-01;
         2.1740265817E-01;2.1945184443E-01;2.2355647425E-01;2.2973339499E-01;2.4200234966E-01;
         2.4409812160E-01;2.7134111043E-01;2.9250000000E-01;]/H3(1);

tw_d3 = tw_d3 - tw_d3(1);


Epl_d3 = [-1.1427949963E+01;-1.1427666484E+01;-1.1427520588E+01;-1.1426767810E+01;-1.1426251000E+01;
          -1.1425637465E+01;-1.1424132754E+01;-1.1421185049E+01;-1.1420022093E+01;-1.1417420572E+01;
          -1.1412838976E+01;-1.1409290655E+01;-1.1407479667E+01;-1.1403458422E+01;-1.1396836126E+01;
          -1.1387035451E+01;-1.1384387193E+01;-1.1378889427E+01;-1.1370159292E+01;-1.1351454291E+01;
          -1.1348148631E+01;-1.1307087178E+01;-1.1343200734E+01;] * Ha2eV * Nc3/H3(1);

rEpl_d3 = Epl_d3 - Epl_d3(1);

ppval(fnder(spline(tw_d3,rEpl_d3),2),tw_d3);%12

% (16,0) CNT
%----------------

% Radius = 11.77 Bohr, mesh-size=0.15, N_c = 16, N_eta = 23, H = 4.0007020478 Bohr

H4 = [4.0009279109;4.0010575997;4.0018321347;4.0040964347;4.0049021242;
      4.0069102595;4.0108751750;4.0124205113;4.0163449648;4.0283632334;
      4.0314105691;4.0372166619;4.0480815878;4.0628574309;4.0662746618;
      4.0749561671;4.0902459493;4.1230483591;4.1284582813;4.1407489517;
      4.2253852808;4.3534365495;4.8072168842;]/nm2Bohr;

st4 = (H4 - H4(1))*100/H4(1);


Nc4 = 16;
ndata4 = 12;

tw_d4 = [1.9634954085E-01;2.0035207865E-01;2.0235292756E-01;2.0835643072E-01;2.1035878003E-01;
         2.1436424064E-01;2.2037084077E-01;2.2237634281E-01;2.2638887772E-01;2.3641341661E-01;
         2.3842736415E-01;2.4245779955E-01;2.4851177672E-01;2.5655138854E-01;2.5858173984E-01;
         2.6264652848E-01;2.6875634020E-01;2.8092166843E-01;2.8298298715E-01;2.8711141209E-01;
         3.0977822516E-01;3.4358130741E-01;3.7260000000E-01;]/H4(1);

tw_d4 = tw_d4 - tw_d4(1);


Epl_d4 = [-1.1426711309E+01;-1.1426598230E+01;-1.1426448438E+01;-1.1425629264E+01;-1.1425232103E+01;
          -1.1424254685E+01;-1.1422340962E+01;-1.1421580175E+01;-1.1419879437E+01;-1.1414619567E+01;
          -1.1413395081E+01;-1.1410781294E+01;-1.1406463724E+01;-1.1400012573E+01;-1.1398263346E+01;
          -1.1394609821E+01;-1.1388784569E+01;-1.1376010672E+01;-1.1373704893E+01;-1.1368974610E+01;
          -1.1340432724E+01;-1.1300879307E+01;-1.1342577063E+01;]* Ha2eV * Nc4/H4(1); % input is in eV


rEpl_d4 = Epl_d4 - Epl_d4(1);

ppval(fnder(spline(tw_d4,rEpl_d4),2),tw_d4);

% (27,0) CNT
%----------------

% Radius = 19.84 Bohr, mesh-size=0.15, N_c = 27, N_eta = 23, H = 4.0025688504 Bohr

H5 = [4.0009179631;4.0030093866;4.0037325087;4.0065547545;4.0096945517;
      4.0117115723;4.0177753721;4.0276339323;4.0333589895;4.0433888941;
      4.0639158055;4.0769812762;4.0854462633;4.1024528857;4.1323366745;
      4.1752159167;4.1874854861;4.2226060302;4.2699963049;4.3496045698;
      4.7958390581;]/nm2Bohr;

st5 = (H5 - H5(1))*100/H5(1);

Nc5 = 27;
ndata5 = 9;

tw_d5 = [1.1635528347E-01;1.2035241978E-01;1.2235392043E-01;1.2635773064E-01;1.2834874714E-01;
         1.3035365326E-01;1.3436573243E-01;1.4036067920E-01;1.4237433968E-01;1.4640616419E-01;
         1.5246922574E-01;1.5647845252E-01;1.5851678719E-01;1.6260121166E-01;1.6875485897E-01;
         1.7686626873E-01;1.7895316195E-01;1.8314101944E-01;1.8947237551E-01;2.0275574747E-01;
         2.1600000000E-01;]/H5(1);

tw_d5 = tw_d5 - tw_d5(1);


Epl_d5 = [-1.1428967467E+01;-1.1428662783E+01;-1.1428234183E+01;-1.1426860899E+01;-1.1425873151E+01;
          -1.1424806978E+01;-1.1422076879E+01;-1.1416762552E+01;-1.1414657124E+01;-1.1409958378E+01;
          -1.1401766975E+01;-1.1395641466E+01;-1.1392334535E+01;-1.1385333823E+01;-1.1373947016E+01;
          -1.1357525086E+01;-1.1353088612E+01;-1.1344207916E+01;-1.1331128288E+01;-1.1304571139E+01;
          -1.1342271743E+01;] * Ha2eV * Nc5/H5(1);


rEpl_d5 = Epl_d5 - Epl_d5(1);

ppval(fnder(spline(tw_d5,rEpl_d5),2),tw_d5);%5

% (22,11) CNT
%-------------

% Radius = 21.4 Bohr, mesh-size=0.15, N_c = 11, N_eta = 30, H = 1.51 Bohr

H6 = [1.5118;1.512;1.5123633431;1.513;1.513;
      1.515;1.5154242765;1.5176956813;1.5200128682;1.5237162234;
      1.5235728247;1.5367542166;1.5398988960;1.5491756480;1.5633193401;
      1.5750998540;1.5890922223;1.6555675223;2.8504144190;]/nm2Bohr;

st6 = (H6 - H6(1))*100/H6(1);

Nc6 = 11;
ndata6 = 11;

tw_d6 = [2.0399952296E-01;2.0490711765E-01;2.0566365359E-01;2.0672337597E-01;2.0748026143E-01;
         2.0899438775E-01;2.0944901399E-01;2.1020652715E-01;2.1172253688E-01;2.1309469021E-01;
         2.1385621106E-01;2.1767549311E-01;2.1844519508E-01;2.1998861680E-01;2.2231741249E-01;
         2.2402617734E-01;2.2561234863E-01;2.2987851440E-01;2.3400000000E-01;]/H6(1);

tw_d6 = tw_d6 - tw_d6(1);

Epl_d6 = [-1.1429039619E+01;-1.1428983091E+01;-1.1428585673E+01;-1.1427715694E+01;-1.1427013974E+01;
          -1.1424581521E+01;-1.1423740257E+01;-1.1422266790E+01;-1.1418503001E+01;-1.1414404128E+01;
          -1.1411824463E+01;-1.1395782700E+01;-1.1392025001E+01;-1.1383893520E+01;-1.1370304963E+01;
          -1.1359467809E+01;-1.1348843615E+01;-1.1315844367E+01;-1.1344022252E+01;]  * Ha2eV * Nc6/H6(1);

rEpl_d6 = Epl_d6 - Epl_d6(1);

ppval(fnder(spline(tw_d6,rEpl_d6),2),tw_d6);%4


% (18,9) CNT
%-----------

% Radius = 17.46 Bohr, mesh-size=0.15, N_c = 9, N_eta = 30, H = 1.51 Bohr

H7 = [1.512;1.5125;1.5125;1.513;1.513;
      1.5134083461;1.514;1.5144881234;1.5161607846;1.5182694689;
      1.5193912835;1.5235632459;1.5274453801;1.5292837116;1.5340639341;
      1.5424144510;1.5460597835;1.5500654725;1.5571407841;1.5710354972;
      1.5837031685;1.5897466373;1.6030806613;1.6838621551;3.1534969490;]/nm2Bohr;

st7 = (H7 - H7(1))*100/H7(1);

Nc7 = 9;
ndata7 = 12;

tw_d7 = [2.4933275028E-01;2.5024084061E-01;2.5099772575E-01;2.5205721665E-01;2.5281415238E-01;
         2.5432826203E-01;2.5478274034E-01;2.5554042055E-01;2.5705642560E-01;2.5842160018E-01;
         2.5918163688E-01;2.6070293042E-01;2.6298522948E-01;2.6375043228E-01;2.6528305966E-01;
         2.6758942451E-01;2.6850166271E-01;2.6927766352E-01;2.7083382340E-01;2.7318125042E-01;
         2.7503159466E-01;2.7582903300E-01;2.7743044938E-01;2.8270618938E-01;2.8800000000E-01;]/H7(1);

tw_d7 = tw_d7 - tw_d7(1);

Epl_d7 = [-1.1428525490E+01;-1.1428373053E+01;-1.1428176167E+01;-1.1427619198E+01;-1.1427126853E+01;
          -1.1425607775E+01;-1.1424980632E+01;-1.1423993358E+01;-1.1421447083E+01;-1.1418688763E+01;
          -1.1416959141E+01;-1.1413071827E+01;-1.1406249692E+01;-1.1403700614E+01;-1.1398174606E+01;
          -1.1388922336E+01;-1.1384990610E+01;-1.1381458798E+01;-1.1374062816E+01;-1.1362194020E+01;
          -1.1352240013E+01;-1.1347829553E+01;-1.1338747292E+01;-1.1305733761E+01;-1.1343591165E+01;]  * Ha2eV * Nc7/H7(1);

rEpl_d7 = Epl_d7 - Epl_d7(1);

ppval(fnder(spline(tw_d7,rEpl_d7),2),tw_d7);%7

% Twist energy Vs twist plot
%------------

figure1 = figure;
%ax=gca;
%ax.Position = [0.13 0.13 1.55 0.8];
box on

hold on
%plot(tw_d2,rEpl_d2,'Marker','o','Markersize',10,'LineWidth',1,'LineStyle','-')

plot(tw_d4,rEpl_d4,'color',[0.4940, 0.1840, 0.5560],'Marker','o','Markersize',10,'LineWidth',1,'LineStyle','-');
plot(tw_d3,rEpl_d3,'color',[0 0.6 0],'Marker','p','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d1,rEpl_d1,'color',[0 0 1],'Marker','d','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d7,rEpl_d7,'color',[1 0 0],'Marker','h','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d5,rEpl_d5,'color',[0.64 0.08 0.18],'Marker','sq','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d6,rEpl_d6,'color',[0.95 0.00 0.75],'Marker','^','Markersize',10,'LineWidth',1,'LineStyle','-')


% S1=scatter(tw_d4(4),rEpl_d4(4),144,'d','Marker','o');S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560]; S1.MarkerFaceColor = [0.4940, 0.1840, 0.5560];
% S2=scatter(tw_d3(4),rEpl_d3(4),144,'d','Marker','p');S2.MarkerEdgeColor = [0 0.6 0]; S2.MarkerFaceColor = [0 0.6 0];
% S3=scatter(tw_d1(4),rEpl_d1(4),144,'d','Marker','d');S3.MarkerEdgeColor = [0 0 1]; S3.MarkerFaceColor = [0 0 1];
% S4=scatter(tw_d7(4),rEpl_d7(4),144,'d','Marker','h');S4.MarkerEdgeColor = [1 0 0]; S4.MarkerFaceColor = [1 0 0];
% S5=scatter(tw_d5(4),rEpl_d5(4),144,'d','Marker','sq');S5.MarkerEdgeColor = [0.64 0.08 0.18]; S5.MarkerFaceColor = [0.64 0.08 0.18];
% S6=scatter(tw_d6(4),rEpl_d6(4),144,'d','Marker','^');S6.MarkerEdgeColor = [0.95 0.00 0.75]; S6.MarkerFaceColor = [0.95 0.00 0.75];



set(gca,'TickLabelInterpreter','LaTex');
set(gca,'yscale','log','XTick',[0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80],'XMinorTick','on','YTick',[1e-1 1e0 1e1 1e2 1e3],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
%set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
set(gca,'YTickLabel',{'$10^{-1}$', '$10^{0}$', '$10^{1}$','$10^{2}$','$10^{3}$'});
xlim([0.00 0.85]);
ylim([1e-1 1.3e3]);

xlabel('Applied twist (rad/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Strain energy density (eV/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');




yyaxis right

plot(tw_d4(5:end),st4(5:end),'color',[0.4940, 0.1840, 0.5560],'Marker','o','MarkerFaceColor',[0.4940, 0.1840, 0.5560],'Markersize',10,'LineWidth',1,'LineStyle','--')
plot(tw_d3(3:end),st3(3:end),'color',[0 0.6 0],'Marker','p','MarkerFaceColor',[0 0.6 0],'Markersize',10,'LineWidth',1,'LineStyle','--')
plot(tw_d1(5:end),st1(5:end),'color',[0 0 1],'Marker','d','MarkerFaceColor',[0 0 1],'Markersize',10,'LineWidth',1,'LineStyle','--')
plot(tw_d7(6:end),st7(6:end),'color',[1 0 0],'Marker','h','MarkerFaceColor',[1 0 0],'Markersize',10,'LineWidth',1,'LineStyle','--')
plot(tw_d5(3:end),st5(3:end),'color',[0.64, 0.08, 0.18],'Marker','sq','MarkerFaceColor',[0.64, 0.08, 0.18],'Markersize',10,'LineWidth',1,'LineStyle','--')
plot(tw_d6(5:end),st6(5:end),'color',[0.95 0.00 0.75],'Marker','^','MarkerFaceColor',[0.95 0.00 0.75],'Markersize',10,'LineWidth',1,'LineStyle','--')


set(gca,'TickLabelInterpreter','LaTex');
set(gca,'yscale','log','Ycolor',[0 0 0],'YTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2],'XTick',[0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80],'XMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
%set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
%set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
set(gca,'YTickLabel',{'$10^{-3}$','$10^{-2}$', '$10^{-1}$', '$10^{0}$','$10^{1}$','$10^{2}$'});

ylim([6e-2 120]);

ylabel('Axial strain (\%)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');


%legend1 = legend('$d_t = 1.25$ nm','$d_t = 1.56$ nm','$d_t = 1.75$ nm','$d_t = 1.85$ nm','$d_t = 2.10$ nm','$d_t = 2.26$ nm','Location','Southeast');
legend1 = legend('$(16,0)$','$(20,0)$','$(13,13)$','$(18,9)$','$(27,0)$','$(22,11)$','Orientation','horizontal');
set(legend1,'Position',[0.05 0.15 1 0.05],'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');



set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPosition',[pos(1) pos(2) pos(3)*2 pos(4)*1.2],'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'strainenergy','epsc')
hold off

% figure1_inset = twistenergyVstwist_linear_relaxed;

% fig3=inset(figure1,figure1_inset);
%  saveas(gca,'twistplot','epsc')

% Strain with twist
figure2 = figure;
box on

hold on
%plot(tw_d2,st2,'Marker','o','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d4(6:end),st4(6:end),'color',[0.4940, 0.1840, 0.5560],'Marker','o','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d3(3:end),st3(3:end),'color',[0 0.6 0],'Marker','p','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d1(6:end),st1(6:end),'color',[0 0 1],'Marker','d','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d7(7:end),st7(7:end),'color',[1 0 0],'Marker','h','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d5(4:end),st5(4:end),'color',[0.64, 0.08, 0.18],'Marker','sq','Markersize',10,'LineWidth',1,'LineStyle','-')
plot(tw_d6(6:end),st6(6:end),'color',[0.95 0.00 0.75],'Marker','^','Markersize',10,'LineWidth',1,'LineStyle','-')

% S1=scatter(tw_d4(4),st4(4),144,'d','Marker','o');S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560]; S1.MarkerFaceColor = [0.4940, 0.1840, 0.5560];
% S2=scatter(tw_d3(4),st3(4),144,'d','Marker','p');S2.MarkerEdgeColor = [0 0.6 0]; S2.MarkerFaceColor = [0 0.6 0];
% S3=scatter(tw_d1(4),st1(4),144,'d','Marker','d');S3.MarkerEdgeColor = [0 0 1]; S3.MarkerFaceColor = [0 0 1];
% S4=scatter(tw_d7(4),st7(4),144,'d','Marker','h');S4.MarkerEdgeColor = [1 0 0]; S4.MarkerFaceColor = [1 0 0];
% S5=scatter(tw_d5(4),st5(4),144,'d','Marker','sq');S5.MarkerEdgeColor = [0.64 0.08 0.18]; S5.MarkerFaceColor = [0.64 0.08 0.18];
% S6=scatter(tw_d6(4),st6(4),144,'d','Marker','^');S6.MarkerEdgeColor = [0.95 0.00 0.75]; S6.MarkerFaceColor = [0.95 0.00 0.75];

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'yscale','log','YTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2],'XTick',[0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80],'XMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
%set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
%set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
set(gca,'YTickLabel',{'$10^{-3}$','$10^{-2}$', '$10^{-1}$', '$10^{0}$','$10^{1}$','$10^{2}$'});

ylim([1e-1 115]);
xlim([0.00 0.85]);

ylabel('Axial strain (\%)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
xlabel('Applied twist (rad/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

%legend2 = legend('$d_t = 1.25$ nm','$d_t = 1.56$ nm','$d_t = 1.75$ nm','$d_t = 1.85$ nm','$d_t = 2.10$ nm','$d_t = 2.26$ nm','Location','Southeast');
%legend2 = legend('$(16,0)$','$(20,0)$','$(13,13)$','$(18,9)$','$(27,0)$','$(22,11)$','Location','South','Orientation','horizontal');
%set(legend2,'Position',[0.05 0.2 1 0.05],'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');


% axes('parent',figure2,'position',[0.30 0.25 0.3 0.25]);
% box on
% hold on
% plot(tw_d4(1:5),st4(1:5),'color',[0.4940, 0.1840, 0.5560],'Marker','o','Markersize',5,'LineWidth',1,'LineStyle','-');
% plot(tw_d3(1:5),st3(1:5),'color',[0 0.6 0],'Marker','p','Markersize',5,'LineWidth',1,'LineStyle','-')
% plot(tw_d1(1:5),st1(1:5),'color',[0 0 1],'Marker','d','Markersize',5,'LineWidth',1,'LineStyle','-')
% plot(tw_d7(1:5),st7(1:5),'color',[1 0 0],'Marker','h','Markersize',5,'LineWidth',1,'LineStyle','-')
% plot(tw_d5(1:5),st5(1:5),'color',[0.64 0.08 0.18],'Marker','sq','Markersize',5,'LineWidth',1,'LineStyle','-')
% plot(tw_d6(1:5),st6(1:5),'color',[0.95 0.00 0.75],'Marker','^','Markersize',5,'LineWidth',1,'LineStyle','-')


% set(gca,'TickLabelInterpreter','LaTex');
% set(gca,'yscale','log','XTick',[0.01 0.03 0.05 0.07],'XMinorTick','on','YTick',[1e-3 1e-2 1e-1],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',12);
% set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
% %set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
% set(gca,'YTickLabel',{'$10^{-3}$', '$10^{-2}$', '$10^{-1}$'});
% xlim([0.01 0.07]);
% ylim([1e-3 0.15]);



%str = 'Slope = 2.0';
%dim = [.2 .58 .3 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

set(figure2,'Units','Inches');
pos2 = get(figure2,'Position');
set(figure2,'PaperPosition',[pos2(1) pos2(2) pos2(3)*2 pos2(4)*0.8],'PaperPositionMode','manual','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'poyntingeffect','epsc')
hold off


%%% With and without relaxation
% (27,0) CNT (w/o relaxation)
%----------------

% H5_nr = 4.0019793139/nm2Bohr;
% Nc5_nr = 27;
% ndata5_nr = 26;
% rate5_nr = 0.0004;
% tw_d5_nr = zeros(ndata5_nr,1);
% for i = 2:ndata5_nr
%     tw_d5_nr(i) = tw_d5_nr(i-1) + (i-1)*rate5_nr * nm2Bohr; % input is in rad/nm
% end

% Epl_d5_nr = [ 
%        -1.1429003617E+01
%        -1.1428945290E+01
%        -1.1428471244E+01
%        -1.1426863582E+01
%        -1.1423073215E+01
%        -1.1415723670E+01
%        -1.1403019472E+01
%        -1.1382970700E+01
%        -1.1353554459E+01
%        -1.1312772585E+01
%        -1.1258115348E+01
%        -1.1189385725E+01
%        -1.1109112093E+01
%        -1.1015257343E+01
%        -1.0913889889E+01
%        -1.0823461361E+01
%        -1.0704047494E+01
%        -1.0527637779E+01
%        -1.0310597768E+01
%        -1.0101068447E+01
%        -9.9734415898E+00
%        -9.9944201702E+00
%        -1.0173976076E+01
%        -1.0448728114E+01
%        -1.0710775185E+01
%        -1.0876529873E+01 
% ] * Ha2eV * Nc5_nr/H5_nr; % input is in eV/nm

% rEpl_d5_nr = Epl_d5_nr - Epl_d5_nr(1);

% figure3 = figure;
% box on

% hold on
% plot(tw_d5,rEpl_d5,'color',[1 0 0],'Marker','sq','Markersize',12,'LineWidth',1,'LineStyle','-')
% plot(tw_d5_nr,rEpl_d5_nr,'color',[0 0 1],'Marker','d','Markersize',12,'LineWidth',1,'LineStyle','-')


% set(gca,'TickLabelInterpreter','LaTex');
% set(gca,'yscale','log','XTick',[0.00 0.50 1.00 1.50 2.00 2.50],'XMinorTick','on','YTick',[1e-1 1e0 1e1 1e2 1e3],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
% set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
% %set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
% %set(gca,'XTickLabel',{'$0.01$', '$0.02$', '$0.04$'});
% %set(gca,'YTickLabel',{'$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
% xlim([0.00 2.50]);
% ylim([1e-1 6e3]);

% xlabel('Twist (rad/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
% ylabel('$\mathcal{E}_{twist}$ (eV/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

% %legend1 = legend('$d_t = 1.25$ nm','$d_t = 1.56$ nm','$d_t = 1.75$ nm','$d_t = 1.85$ nm','$d_t = 2.10$ nm','$d_t = 2.26$ nm','Location','Southeast');
% legend3 = legend('relaxed','unrelaxed','Location','Southeast');

% set(legend3,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

% %str = 'Slope = 2.0';
% %dim = [.2 .58 .3 .3];
% %annotation('textbox',dim,'String',str,'FitBoxToText','on');

% set(figure3,'Units','Inches');
% pos = get(figure3,'Position');
% set(figure3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
% saveas(gcf,'twistenergyVstwist_comparison','epsc')
% hold off