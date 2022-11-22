close all
clear
clc

nm2Bohr = 18.89726125;
Ha2eV = 27.211386245988;
%indx_ln = 5;

% (13,13) CNT
%----------------

% Radius = 16.56 Bohr, mesh-size=0.125, N_c = 13, N_eta = 40, H = 2.309 Bohr

%H1 = [2.3106970869;2.3092161903;2.3103942073;2.3102709652;2.3105578117;2.3109540520]; %2.3091577619;
H1 = [2.3097769325;2.3107113331;2.3099969252;2.3102540474;2.3105781717;2.3102130272]/nm2Bohr;
Nc1 = 13;
ndata1 = 5;
rate1 = 0.0004;
%tw_d1= [2.4166097335E-01; 2.4258525219E-01;2.4350893866E-01;2.4443309635E-01;2.4535720473E-01;2.4628142786E-01]; % input is in rad i.e. total twist per unit cell
tw_d1= [2.4166097335E-01
        %2.4258488413E-01
        2.4350916866E-01
        2.4443316743E-01
        2.4535726905E-01
        2.4628150032E-01]/H1(1);

tw_d1_f= [2.4166097335E-01
        %2.4258488413E-01
        2.4350916866E-01
        2.4443316743E-01
        2.4535726905E-01
        2.4628150032E-01]/H1(1);        

% energy at 0.15 mesh
% Epl_d1 = [
%           -1.1428378403E+01
%           -1.1428384137E+01
%           -1.1428253882E+01
%           -1.1428075149E+01
%           -1.1427808039E+01
%           -1.1427461148E+01
%  ] * 2 * Ha2eV * Nc1; % input is in Ha -> eV

% energy at 0.125 mesh
% Epl_d1 = [
%           -1.1428407560E+01
%           -1.1428363991E+01
%           -1.1428257194E+01
%           -1.1428068379E+01
%           -1.1427803087E+01
%           -1.1427463566E+01
%  ] * 2 * Ha2eV * Nc1; % input is in Ha -> eV

% energy at 0.1 mesh
Epl_d1 = [
          -1.1428411403E+01
          %-1.1428353001E+01
          -1.1428261184E+01
          -1.1428071709E+01
          -1.1427804755E+01
          -1.1427466760E+01
 ] * 2 * Ha2eV * Nc1/H1(1); % input is in Ha -> eV

 Epl_d1_f = [
          -1.1428411403E+01
          %-1.1428353001E+01
          -1.1428261184E+01
          -1.1428071709E+01
          -1.1427804755E+01
          -1.1427466760E+01
 ] * 2 * Ha2eV * Nc1/H1(1);

tw_d1 = tw_d1 - tw_d1(1);
tw_d1 = tw_d1(2:end);
rEpl_d1 = Epl_d1 - Epl_d1(1);
rEpl_d1 = rEpl_d1(2:end);

tw_d1_f = tw_d1_f - tw_d1_f(1);
tw_d1_f = tw_d1_f(2:end);
tw2_d1_f = tw_d1_f.^2;
rEpl_d1_f = Epl_d1_f - Epl_d1_f(1);
rEpl_d1_f = rEpl_d1_f(2:end);

log_tw_d1 = log10(tw_d1);%(1:indx_ln)
log_rEpl_d1 = log10(rEpl_d1);
P_Efit_d1 = polyfit(log_tw_d1,log_rEpl_d1,1);
Efit_d1 = P_Efit_d1(1)*log_tw_d1 + P_Efit_d1(2);
antilog_Efit_d1= 10.^(Efit_d1);


% (15,0) CNT
%----------------

% Radius = 11 Bohr, mesh-size=0.125, N_c = 15, N_eta = 23, H = 4.0017685908 Bohr

%H2 = [4.0017691813;4.0018078329;4.0019438369;4.0021520915;4.0022985527;4.0027004324];
H2 = [4.0018117784;4.0018514174;4.0019626788;4.0021865072;4.0023588875;4.0027583226]/nm2Bohr;
Nc2 = 15;
ndata2 = 5;
rate2 = 0.0004;
%tw_d2= [2.0943951024E-01;2.1104021791E-01;2.1264094104E-01;2.1424171858E-01;2.1584257942E-01;2.1744349884E-01]; % input is in rad
tw_d2 = [2.0943951024E-01
         %2.1104023495E-01
         2.1264097552E-01
         2.1424176059E-01
         2.1584263519E-01
         2.1744357875E-01]/H2(1);
         

% energy at 0.15 mesh
% Epl_d2 = [
%           -1.1426118796E+01
%           -1.1426102722E+01
%           -1.1426051529E+01
%           -1.1425965304E+01
%           -1.1425844101E+01
%           -1.1425688211E+01
% ] * 2 * Ha2eV * Nc2 ; % input is in Ha -> eV

% % energy at 0.125 mesh
% Epl_d2 = [
%           -1.1426085685E+01
%           -1.1426070620E+01
%           -1.1426021042E+01
%           -1.1425936628E+01
%           -1.1425818835E+01
%           -1.1425665114E+01
% ] * 2 * Ha2eV * Nc2 ; % input is in Ha -> eV

% energy at 0.1 mesh
Epl_d2 = [
          -1.1426088918E+01
          %-1.1426074127E+01
          -1.1426024755E+01
          -1.1425940407E+01
          -1.1425822138E+01
          -1.1425668261E+01
] * 2 * Ha2eV * Nc2/H2(1) ; % input is in Ha -> eV

rEpl_d2 = Epl_d2 - Epl_d2(1);
tw_d2 = tw_d2 - tw_d2(1);
tw_d2 = tw_d2(2:end);
tw2_d2 = tw_d2.^2;
rEpl_d2 = rEpl_d2(2:end);

log_tw_d2 = log10(tw_d2);
log_rEpl_d2 = log10(rEpl_d2);
P_Efit_d2 = polyfit(log_tw_d2,log_rEpl_d2,1);
Efit_d2 = P_Efit_d2(1)*log_tw_d2 + P_Efit_d2(2);
antilog_Efit_d2= 10.^(Efit_d2);



% (20,0) CNT
%----------------

% Radius = 14.72 Bohr, mesh-size=0.125, N_c = 20, N_eta = 23, H = 4.0002688791 Bohr

%H3 = [4.0002709862;4.0004777268;4.0007616136;4.0006818012;4.0014526732;4.0019124672];
H3 = [4.0012612687;4.0002692946;4.0006929657;4.0026067611;4.0015353019;4.0019658895]/nm2Bohr;
Nc3 = 20;
ndata3 = 5;
rate3 = 0.0004;
%tw_d3 = [1.5707963268E-01;1.5867974107E-01;1.6027993216E-01;1.6188023681E-01;1.6348050953E-01;1.6508109060E-01]; % input is in rad 
tw_d3 = [1.5707963268E-01
         %1.5868013719E-01
         1.6028024490E-01
         1.6188052209E-01
         1.6348156480E-01
         1.6508217892E-01]/H3(1);

tw_d3_f = [1.5707963268E-01
         %1.5868013719E-01
         1.6028024490E-01
         1.6188052209E-01
         1.6348156480E-01
         1.6508217892E-01]/H3(1);         

% energy at 0.15 mesh
% Epl_d3 = [-1.1427964269E+01
%           %-1.1427932827E+01
%           -1.1427839042E+01
%           -1.1427678005E+01
%           -1.1427465542E+01
%           -1.1427186913E+01
% ] * 2 * Ha2eV * Nc3 ; % input is in Ha -> eV

% energy at 0.125 mesh
% Epl_d3 = [
%           -1.1427932649E+01
%           -1.1427905489E+01
%           -1.1427811956E+01
%           -1.1427651445E+01
%           -1.1427444875E+01
%           -1.1427170993E+01
% ] * 2 * Ha2eV * Nc3 ; % input is in Ha -> eV

% energy at 0.1 mesh
Epl_d3 = [
         -1.1427937009E+01
         %-1.1427907076E+01
         -1.1427815658E+01
         -1.1427654463E+01
         -1.1427448996E+01
         -1.1427174220E+01
] * 2 * Ha2eV * Nc3/H3(1) ; % input is in Ha -> eV

Epl_d3_f = [
         -1.1427937009E+01
         %-1.1427907076E+01
         -1.1427815658E+01
         -1.1427654463E+01
         -1.1427448996E+01
         -1.1427174220E+01
] * 2 * Ha2eV * Nc3/H3(1) ;

rEpl_d3 = Epl_d3 - Epl_d3(1);
tw_d3 = tw_d3 - tw_d3(1);
tw_d3 = tw_d3(2:end);
rEpl_d3 = rEpl_d3(2:end);


tw_d3_f = tw_d3_f - tw_d3_f(1);
tw_d3_f = tw_d3_f(2:end);
tw2_d3_f = tw_d3_f.^2;
rEpl_d3_f = Epl_d3_f - Epl_d3_f(1);
rEpl_d3_f = rEpl_d3_f(2:end);


log_tw_d3 = log10(tw_d3);
log_rEpl_d3 = log10(rEpl_d3);
P_Efit_d3 = polyfit(log_tw_d3,log_rEpl_d3,1);
Efit_d3 = P_Efit_d3(1)*log_tw_d3 + P_Efit_d3(2);
antilog_Efit_d3= 10.^(Efit_d3);


% (16,0) CNT
%----------------

% Radius = 11.78 Bohr, mesh-size=0.15, N_c = 16, N_eta = 23, H = 4.0023317752 Bohr

H4 = [4.0023317752;4.0021629147;4.0027071210E+00;4.0026260828;4.0029036314;4.0032584275]/nm2Bohr;
Nc4 = 16;
ndata4 = 4;
rate4 = 0.0004;
tw_d4 = [1.9634954085E-01;2.0115242157E-01;2.0275347201E-01;2.0435463346E-01]/H4(1); % input is in rad;
%1.9795047356E-01;1.9955133873E-01;

tw_d4_f = [1.9634954085E-01;1.9955133873E-01;2.0115242157E-01;2.0275347201E-01;2.0435463346E-01]/H4(1); % input is in rad;
% energy at 0.15 mesh
% Epl_d4 = [
%           -1.1426612656E+01
%           -1.1426592801E+01
%           -1.1426534853E+01
%           -1.1426438287E+01
%           -1.1426302562E+01
%           -1.1426127969E+01
% ] * 2 * Ha2eV * Nc4 ; % input is in Ha -> eV

% energy at 0.1 mesh
Epl_d4 = [
          -1.1426624938E+01
          %-1.1426607368E+01
          %-1.1426546657E+01
          -1.1426452281E+01
          -1.1426316776E+01
          -1.1426142391E+01
] * 2 * Ha2eV * Nc4/H4(1) ; % input is in Ha -> eV

Epl_d4_f = [
          -1.1426624938E+01
          %-1.1426607368E+01
          -1.1426546657E+01
          -1.1426452281E+01
          -1.1426316776E+01
          -1.1426142391E+01
] * 2 * Ha2eV * Nc4/H4(1) ; % input is in Ha -> eV

rEpl_d4 = Epl_d4 - Epl_d4(1);
tw_d4 = tw_d4 - tw_d4(1);
tw_d4 = tw_d4(2:end);
tw2_d4 = tw_d4.^2;
rEpl_d4 = rEpl_d4(2:end);

tw_d4_f = tw_d4_f - tw_d4_f(1);
tw_d4_f = tw_d4_f(2:end);
tw2_d4_f = tw_d4_f.^2;
rEpl_d4_f = Epl_d4_f - Epl_d4_f(1);
rEpl_d4_f = rEpl_d4_f(2:end);

log_tw_d4 = log10(tw_d4);
log_rEpl_d4 = log10(rEpl_d4);
P_Efit_d4 = polyfit(log_tw_d4,log_rEpl_d4,1);
Efit_d4 = P_Efit_d4(1)*log_tw_d4 + P_Efit_d4(2);
antilog_Efit_d4= 10.^(Efit_d4);


% (27,0) CNT
%----------------

% Radius = 19.84 Bohr, mesh-size=0.15, N_c = 27, N_eta = 23, H = 4.0019793139 Bohr

H5 = [4.0019793139;4.0041416052;4.0030856480;4.0043754833;4.0062230584;4.0050847336]/nm2Bohr;
Nc5 = 27;
ndata5 = 5;
rate5 = 0.0004;
tw_d5 = [1.1635528347E-01;1.1955773183E-01;1.2115896609E-01;1.2276071629E-01;1.2436320551E-01]/H5(1); % input is in rad;
%1.1795607519E-01;
tw_d5_f = [1.1635528347E-01;1.1955773183E-01;1.2115896609E-01;1.2276071629E-01;1.2436320551E-01]/H5(1); % input is in rad;
% energy at 0.15 mesh
% Epl_d5 = [
%           -1.1428970228E+01
%           -1.1428858759E+01
%           -1.1428730206E+01
%           -1.1428430356E+01
%           -1.1428019445E+01
%           -1.1427573325E+01
% ] * 2 * Ha2eV * Nc5 ; % input is in Ha -> eV

% energy at 0.1 mesh
Epl_d5 = [
          -1.1428997073E+01
          %-1.1428911774E+01
          -1.1428771314E+01
          -1.1428483236E+01
          -1.1428072848E+01
          -1.1427607354E+01         
] * 2 * Ha2eV * Nc5/H5(1) ; % input is in Ha -> eV

Epl_d5_f = [
          -1.1428997073E+01
          %-1.1428911774E+01
          -1.1428771314E+01
          -1.1428483236E+01
          -1.1428072848E+01
          -1.1427607354E+01         
] * 2 * Ha2eV * Nc5/H5(1) ; % input is in Ha -> eV


rEpl_d5 = Epl_d5 - Epl_d5(1);
tw_d5 = tw_d5 - tw_d5(1);
tw_d5 = tw_d5(2:end);
tw2_d5 = tw_d5.^2;
rEpl_d5 = rEpl_d5(2:end);

tw_d5_f = tw_d5_f - tw_d5_f(1);
tw_d5_f = tw_d5_f(2:end);
tw2_d5_f = tw_d5_f.^2;
rEpl_d5_f = Epl_d5_f - Epl_d5_f(1);
rEpl_d5_f = rEpl_d5_f(2:end);

log_tw_d5 = log10(tw_d5);
log_rEpl_d5 = log10(rEpl_d5);
P_Efit_d5 = polyfit(log_tw_d5,log_rEpl_d5,1);
Efit_d5 = P_Efit_d5(1)*log_tw_d5 + P_Efit_d5(2);
antilog_Efit_d5= 10.^(Efit_d5);

% (22,11) CNT
%----------------

% Radius = 21.38 Bohr, mesh-size=0.15, N_c = 11, N_eta = 70, H = 1.5 Bohr

H6 = [1.5124912124;1.5125777267;1.5129790920;1.5137618472;1.5135718540;1.5138339053]/nm2Bohr;
Nc6 = 11;
ndata6 = 5;
rate6 = 0.0004;
tw_d6 = [2.0399952296E-01;2.0581497183E-01;2.0642016346E-01;2.0702566820E-01]/H6(1); %;2.0460451945E-01;2.0763109694E-01

tw_d6_f = [2.0399952296E-01;2.0460451945E-01;2.0581497183E-01;2.0642016346E-01;2.0702566820E-01;]/H6(1); %;2.0460451945E-01;2.0763109694E-01
% energy at 0.15 mesh
% Epl_d6 = [
%           -1.1429162599E+01
%           -1.1429092110E+01
%           -1.1428562680E+01
%           -1.1428072827E+01
%           -1.1427524937E+01
%           -1.1426801538E+01
% ] * 2 * Ha2eV * Nc6 ; % input is in Ha -> eV

% energy at 0.1 mesh
Epl_d6 = [
          -1.1429172748E+01
          %-1.1429106044E+01
          -1.1428583391E+01
          -1.1428105894E+01
          -1.1427544348E+01
          %-1.1426829761E+01
 ] * 2 * Ha2eV * Nc6/H6(1) ; % input is in Ha -> eV

 Epl_d6_f = [
          -1.1429172748E+01
          -1.1429106044E+01
          -1.1428583391E+01
          -1.1428105894E+01
          -1.1427544348E+01
          %-1.1426829761E+01
 ] * 2 * Ha2eV * Nc6/H6(1) ; % input is in Ha -> eV

rEpl_d6 = Epl_d6 - Epl_d6(1);
tw_d6 = tw_d6 - tw_d6(1);
tw_d6 = tw_d6(2:end);
tw2_d6 = tw_d6.^2;
rEpl_d6 = rEpl_d6(2:end);

tw_d6_f = tw_d6_f - tw_d6_f(1);
tw_d6_f = tw_d6_f(2:end);
tw2_d6_f = tw_d6_f.^2;
rEpl_d6_f = Epl_d6_f - Epl_d6_f(1);
rEpl_d6_f = rEpl_d6_f(2:end);

log_tw_d6 = log10(tw_d6);
log_rEpl_d6 = log10(rEpl_d6);
P_Efit_d6 = polyfit(log_tw_d6,log_rEpl_d6,1);
Efit_d6 = P_Efit_d6(1)*log_tw_d6 + P_Efit_d6(2);
antilog_Efit_d6= 10.^(Efit_d6);

% (18,9) CNT
%----------------

% Radius = 17.46 Bohr, mesh-size=0.15, N_c = 9, N_eta = 70, H = 1.5 Bohr

H7 = [1.5130233554;1.5123736013;1.5126508325;1.5125876811;1.5130391119;1.5132178566;1.5135765153]/nm2Bohr;%1.5123775448;
Nc7 = 9;
ndata7 = 4;
rate7 = 0.0004;
tw_d7 = [
         2.4933275028E-01
         %2.4993795963E-01
         %2.5054291064E-01
         2.5114786009E-01
         2.5175292042E-01
         2.5235795549E-01
         %2.5308816818E-01
         ]/H7(1);

tw_d7_f = [
         2.4933275028E-01
         %2.4993795963E-01
         2.5054291064E-01
         2.5114786009E-01
         2.5175292042E-01
         2.5235795549E-01
         %2.5308816818E-01
         ]/H7(1);         
% energy at 0.15 mesh
% Epl_d7 = [
%           -1.1428547731E+01
%           %-1.1428555894E+01
%           -1.1428431724E+01
%           -1.1428210544E+01
%           -1.1427925651E+01
%           -1.1427527839E+01
% ] * 2 * Ha2eV * Nc7 ; % input is in Ha -> eV

% energy at 0.1 mesh
Epl_d7 = [ 
          %-1.1428582817E+01
          -1.1428599000E+01
          %-1.1428561804E+01
          %-1.1428433278E+01
          -1.1428214704E+01
          -1.1427918871E+01
          -1.1427525621E+01
          %-1.1426943502E+01     
] * 2 * Ha2eV * Nc7/H7(1) ; % input is in Ha -> eV

Epl_d7_f = [ 
          %-1.1428582817E+01
          -1.1428599000E+01
          %-1.1428561804E+01
          -1.1428433278E+01
          -1.1428214704E+01
          -1.1427918871E+01
          -1.1427525621E+01
          %-1.1426943502E+01     
] * 2 * Ha2eV * Nc7/H7(1) ; % input is in Ha -> eV

rEpl_d7 = Epl_d7 - Epl_d7(1);
tw_d7 = tw_d7 - tw_d7(1);
tw_d7 = tw_d7(2:end);
tw2_d7 = tw_d7.^2;
rEpl_d7 = rEpl_d7(2:end);

tw_d7_f = tw_d7_f - tw_d7_f(1);
tw_d7_f = tw_d7_f(2:end);
tw2_d7_f = tw_d7_f.^2;
rEpl_d7_f = Epl_d7_f - Epl_d7_f(1);
rEpl_d7_f = rEpl_d7_f(2:end);

log_tw_d7 = log10(tw_d7);
log_rEpl_d7 = log10(rEpl_d7);
P_Efit_d7 = polyfit(log_tw_d7,log_rEpl_d7,1);
Efit_d7 = P_Efit_d7(1)*log_tw_d7 + P_Efit_d7(2);
antilog_Efit_d7= 10.^(Efit_d7);

% Twist energy Vs twist plot
%------------

% figure1 = figure;
% box on

% hold on
% %scatter(tw_d2,rEpl_d2,'Marker','o');
% S1=plot(tw2_d4,rEpl_d4*0.5,'d','Marker','o','MarkerSize',12); S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560];
% S2=plot(tw2_d3,rEpl_d3*0.5,'d','Marker','p','MarkerSize',12); S2.MarkerEdgeColor = [0 0.6 0];
% S3=plot(tw2_d1,rEpl_d1*0.5,'d','Marker','d','MarkerSize',12); S3.MarkerEdgeColor = [0 0 1];
% S4=plot(tw2_d7,rEpl_d7*0.5,'d','Marker','h','MarkerSize',12); S4.MarkerEdgeColor = [1 0 0];
% S5=plot(tw2_d5,rEpl_d5*0.5,'d','Marker','sq','MarkerSize',12); S5.MarkerEdgeColor = [0.64 0.08 0.18];
% S6=plot(tw2_d6,rEpl_d6*0.5,'d','Marker','^','MarkerSize',12); S6.MarkerEdgeColor = [0.95 0.00 0.75];

% %plot(tw_d2,antilog_Efit_d2,'LineWidth',1,'LineStyle','-');
% plot(tw2_d4,antilog_Efit_d4*0.5,'color',[0.4940, 0.1840, 0.5560],'LineWidth',1,'LineStyle','-');
% plot(tw2_d3,antilog_Efit_d3*0.5,'color',[0 0.6 0],'LineWidth',1,'LineStyle','-');
% plot(tw2_d1,antilog_Efit_d1*0.5,'color',[0 0 1],'LineWidth',1,'LineStyle','-');
% plot(tw2_d7,antilog_Efit_d7*0.5,'color',[1 0 0],'LineWidth',1,'LineStyle','-');
% plot(tw2_d5,antilog_Efit_d5*0.5,'color',[0.64 0.08 0.18],'LineWidth',1,'LineStyle','-');
% plot(tw2_d6,antilog_Efit_d6*0.5,'color',[0.95 0.00 0.75],'LineWidth',1,'LineStyle','-');

% set(gca,'TickLabelInterpreter','LaTex');
% set(gca,'XTick',[0.0 5.0e-4 1.0e-3 1.5e-3],'YTick',[0 1 2 3 4 5 6],'XMinorTick','on','YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
% set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.1f'));
% %set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
% set(gca,'XTickLabel',{'$0.0$', '$0.5$', '$1.0$','$1.5$'});
% set(gca,'YTickLabel',{'$0$','$1$','$2$','$3$','$4$','$5$','$6$'});
% xlim([0.0000 0.0015]);
% ylim([0 6]);

% xlabel('Twist$^2$ (rad$^2$/nm$^2$)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
% ylabel('$\mathcal{E}_{twist}$ (eV/nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

% legend1 = legend('$(16,0)$','$(20,0)$','$(13,13)$','$(18,9)$','$(27,0)$','$(22,11)$','Location','Northwest');
% set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

% %str = 'Slope = 2.0';
% %dim = [.2 .58 .3 .3];
% %annotation('textbox',dim,'String',str,'FitBoxToText','on');

% set(figure1,'Units','Inches');
% pos = get(figure1,'Position');
% set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4 2]);
% saveas(gcf,'twistenergyVstwist_linear','epsc');
% hold off



% Torsion modulus Vs radius^3
dia3 = [1.25;1.56;1.75;1.85;2.1;2.26].^3;
TM = 10.^[P_Efit_d4(2);P_Efit_d3(2);P_Efit_d1(2);P_Efit_d7(2);P_Efit_d5(2);P_Efit_d6(2)]/1000; % in keV nmP_Efit_d7(2)

figure2 = figure;
box on

hold on
S1=plot(dia3(1),TM(1),'d','Marker','o','MarkerSize',12); S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560];
S2=plot(dia3(2),TM(2),'d','Marker','p','MarkerSize',12); S2.MarkerEdgeColor = [0 0.6 0];
S3=plot(dia3(3),TM(3),'d','Marker','d','MarkerSize',12); S3.MarkerEdgeColor = [0 0 1];
S4=plot(dia3(4),TM(4),'d','Marker','h','MarkerSize',12); S4.MarkerEdgeColor = [1 0 0];
S5=plot(dia3(5),TM(5),'d','Marker','sq','MarkerSize',12); S5.MarkerEdgeColor = [0.64 0.08 0.18];
S6=plot(dia3(6),TM(6),'d','Marker','^','MarkerSize',12); S6.MarkerEdgeColor = [0.95 0.00 0.75];
%plot(dia3,TM,'k','LineWidth',1,'LineStyle','-');

cons = sum(TM(:)./dia3(:))/6 ;
TM_proj = dia3*cons;
plot(dia3,TM_proj,'k','LineWidth',1,'LineStyle','-');

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'XTick',[0 3 6 9 12],'XMinorTick','on','YMinorTick','on','YTick',[0 2 4 6 8],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))

xlim([0 12]);
ylim([0 8.5]);

xlabel('Diameter$^3$ (nm$^3$)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Torsional modulus (keV nm)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

legend2 = legend('$(16,0)$','$(20,0)$','$(13,13)$','$(18,9)$','$(27,0)$','$(22,11)$','Location','Southeast');
set(legend2,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');


axes('parent',figure2,'position',[0.20 0.64 .30 0.27]);
box on

hold on
%scatter(tw_d2,rEpl_d2,'Marker','o');
S1=plot(tw2_d4_f,rEpl_d4_f*0.5,'color',[0.4940, 0.1840, 0.5560],'Marker','o','MarkerSize',8); %S1.MarkerEdgeColor = [0.4940, 0.1840, 0.5560];
S2=plot(tw2_d3_f,rEpl_d3_f*0.5,'color',[0 0.6 0],'Marker','p','MarkerSize',8); %S2.MarkerEdgeColor = [0 0.6 0];
S3=plot(tw2_d1_f,rEpl_d1_f*0.5,'color',[0 0 1],'Marker','d','MarkerSize',8); %S3.MarkerEdgeColor = [0 0 1];
S4=plot(tw2_d7_f,rEpl_d7_f*0.5,'color',[1 0 0],'Marker','h','MarkerSize',8); %S4.MarkerEdgeColor = [1 0 0];
S5=plot(tw2_d5_f,rEpl_d5_f*0.5,'color',[0.64, 0.08, 0.18],'Marker','sq','MarkerSize',8);% S5.MarkerEdgeColor = [0.64 0.08 0.18];
S6=plot(tw2_d6_f,rEpl_d6_f*0.5,'color',[0.95 0.00 0.75],'Marker','^','MarkerSize',8); %S6.MarkerEdgeColor = [0.95 0.00 0.75];


set(gca,'TickLabelInterpreter','LaTex');
set(gca,'YTick',[0 2 4 6],'XMinorTick','on','YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
%set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.4f'));
% %set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
set(gca,'XTickLabel',{'$0$', '$5$', '$10$','$15$'});
set(gca,'YTickLabel',{'$0$','$2$','$4$','$6$'});
xlim([0.0 0.0015]);
ylim([0 6.5]);
% ax = gca;
% ax.XAxis.Exponent = -4;
xlabel('$\alpha_a^2$ (rad$^2$/nm$^2$)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\mathcal{E}_\texttt{twist}$ (eV/nm)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');

set(figure2,'Units','Inches');
pos = get(figure2,'Position');
set(figure2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4 2])
saveas(gcf,'TMvsdia','epsc')
hold off