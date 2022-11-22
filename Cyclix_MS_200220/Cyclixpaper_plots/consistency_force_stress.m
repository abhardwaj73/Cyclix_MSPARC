% Consistency between Energy(Ha/atom)/Displacement and forces(Ha/Bohr) at 0.15 Bohr mesh
%---------------------------------------------------------------------------------------
% Note: Perturbation is given in the direction in which force is checked for consistency
close all
clear
clc

% Energy and force on perturbed atom
natom = 2;
N_c = 9; % cyclic group order
ndata = 5; % on either side of reference
fineness = 0.001;
perturb_start = -0.25;
perturb_step = 0.05;
perturb_end = 0.25;
perturb = [perturb_start:perturb_step:perturb_end]';
fine_pert = [perturb_start:fineness:perturb_end]';

% X-direction

Ha_X = [-5.7083921600;-5.7105409837;-5.7121928923;-5.7133644997;-5.7140689351;-5.7143029315;
        -5.7140802880;-5.7134191953;-5.7123522759;-5.7109088831;-5.7091057939]; % in Ha/atom

F_X = [0.095853405418647;0.076232776205868;0.056504254006336;0.036939733963849;0.017935242905762;3.6373772854E-04;
       -1.7314181107E-02;-3.4291727507E-02;-5.0365536537E-02;-6.5433324086E-02;-7.9315644017E-02;]; % in Ha/Bohr 

% Y-direction

Ha_Y = [-5.7079998739;-5.7102109107;-5.7119611808;-5.7132300658;-5.7140255777;-5.7143029315;
        -5.7140512941;-5.7132758248;-5.7119707291;-5.7101300568;-5.7077552039;]; % in Ha/atom 


F_Y = [9.6367393656E-02;7.9152747166E-02;6.0676455483E-02;4.1239812508E-02;2.1047894682E-02;3.7253420372E-04;
       -0.019826004124958;-0.041305998458969;-0.062988903763609;-0.084264453451195;-0.105201121179974;]; % in Ha/Bohr 

% % Z-direction

Ha_Z = [-5.7046141364;-5.7081775401;-5.7108899263;-5.7128097137;-5.7139378441;-5.7143029315;
        -5.7139528174;-5.7128894954;-5.7111365189;-5.7087274856;-5.7056918679]; % in Ha/atom 

F_Z = [1.5967819322E-01;1.2552029017E-01;9.2390359321E-02;6.0885291617E-02;2.9463457297E-02;2.9403088917E-05;
       -2.8560114251E-02;-5.6691802917E-02;-8.3540670644E-02;-1.0901478200E-01;-1.3350491958E-01]; % in Ha/Bohr


% Energy variation and Energy/Displacement

% X-direction
Ha_var_X = (Ha_X(ndata+1) - Ha_X)*natom;% Force on atom = -variation in total energy/ displacment of the atom
pp_X = spline(perturb,Ha_var_X);
Ha_var_fit_X = ppval(pp_X,fine_pert);
p_der_X = fnder(pp_X,1);
Ha_var_per_disp_X = ppval(p_der_X,fine_pert);

% Y-direction
Ha_var_Y = (Ha_Y(ndata+1) - Ha_Y)*natom;% Force on atom = -variation in total energy/ displacment of the atom
pp_Y = spline(perturb,Ha_var_Y);
Ha_var_fit_Y = ppval(pp_Y,fine_pert);
p_der_Y = fnder(pp_Y,1);
Ha_var_per_disp_Y = ppval(p_der_Y,fine_pert);

% % Z-direction
Ha_var_Z = (Ha_Z(ndata+1) - Ha_Z)*natom;% Force on atom = -variation in total energy/ displacment of the atom
pp_Z = spline(perturb,Ha_var_Z);
Ha_var_fit_Z = ppval(pp_Z,fine_pert);
p_der_Z = fnder(pp_Z,1);
Ha_var_per_disp_Z = ppval(p_der_Z,fine_pert);


% Stress
ndata_S = 7; % on either side of reference
fineness_S = 0.001;
perturb_start_S = -0.07;
perturb_step_S = 0.01;
perturb_end_S = 0.07;
perturb_S = [perturb_start_S:perturb_step_S:perturb_end_S];
fine_pert_S = [perturb_start_S:fineness_S:perturb_end_S]';


% Stress

Ha_S = [-5.7116849160;-5.7124075347;-5.7130059897;-5.7134859616;-5.7138516133;-5.7141063070;-5.7142549675;-5.7143029315;
        -5.7142564463;-5.7141197105;-5.7138949913;-5.7135853256;-5.7131939574;-5.7127245410;-5.7121798132 ]; % Ha/atom

Sig = [-1.4171251516;-1.1888015518;-9.6931323461E-01 ;-7.5869537744E-01 ;-5.5604582082E-01;-3.6223682217E-01;-1.7629815019E-01;7.6427520602E-05;
        1.6633098569E-01;3.2609961012E-01;4.8126889776E-01;6.3192776292E-01;7.7600867847E-01;9.1440373591E-01;1.0468788636];


Ha_var_S = (Ha_S-Ha_S(ndata_S+1) )*natom*N_c;% Stress = variation in total energy of circle/ displacment of the cell
pp_S = spline(perturb_S,Ha_var_S);
Ha_var_fit_S = ppval(pp_S,fine_pert_S);
p_der_S = fnder(pp_S,1);
Ha_var_per_disp_S = ppval(p_der_S,fine_pert_S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure1 = figure;
% box on

% hold on
% S1=plot(perturb,Ha_var_X,'rd','Marker','o','MarkerSize',12);
% S2=plot(perturb,Ha_var_Y,'d','Marker','*','MarkerSize',12); S2.MarkerEdgeColor = [0 0.6 0];
% S3=plot(perturb,Ha_var_Z,'bd','Marker','sq','MarkerSize',12);

% plot(fine_pert,Ha_var_fit_X,'color',[1 0 0],'LineWidth',1,'LineStyle','-')
% plot(fine_pert,Ha_var_fit_Y,'color',[0 1 0],'LineWidth',1,'LineStyle','-')
% plot(fine_pert,Ha_var_fit_Z,'color',[0 0 1],'LineWidth',1,'LineStyle','-')


% xlabel('Displacement (Bohr)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
% ylabel('Energy variation (Ha)','FontSize',16,'FontName','Times New Roman','Interpreter','LaTex');
% legend1 = legend('Force$_X$','Force$_Y$','Force$_Z$');
% set(legend1,'fontsize',16,'FontName','Times New Roman','Interpreter','LaTex');

% set(figure1,'Units','Inches');
% pos = get(figure1,'Position');
% set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
% saveas(gcf,'eggbox','epsc')
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure2 = figure;
box on

hold on
S1=plot(perturb,F_X,'rd','Marker','o','MarkerSize',12);
S2=plot(perturb,F_Y,'d','Marker','*','MarkerSize',12); S2.MarkerEdgeColor = [0 0.6 0];
S3=plot(perturb,F_Z,'bd','Marker','sq','MarkerSize',12);

h(1)=plot(fine_pert,Ha_var_per_disp_X,'color',[1 0 0],'LineWidth',1,'LineStyle','-');
h(2)=plot(fine_pert,Ha_var_per_disp_Y,'color',[0 0.6 0],'LineWidth',1,'LineStyle','-');
h(3)=plot(fine_pert,Ha_var_per_disp_Z,'color',[0 0 1],'LineWidth',1,'LineStyle','-');

set( get( get( h(1), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h(2), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h(3), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set(gca,'XTick',[-0.24 -0.12 0.0 0.12 0.24],'YTick',[-0.15 -0.10 -0.05 0.00 0.05 0.10 0.15],'XMinorTick','on','YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.2f'))
xlim([-0.25 0.25]);
ylim([-0.15 0.17]);

xlabel('Perturbation (Bohr)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Force (Ha/Bohr)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');


yyaxis right
S4=plot(perturb_S,Sig,'d','Marker','d','MarkerSize',12); S4.MarkerEdgeColor = [0.95 0.00 0.75];
plot(fine_pert_S,Ha_var_per_disp_S,'color',[0.95 0.00 0.75],'LineWidth',1,'LineStyle','-')

set(gca,'YTick',[-1.5 -1.0 -0.5 0.0 0.5 1.0],'YMinorTick','on','TickLength',[0.018 0.025],'Ycolor',[0 0 0],'FontName','Times New Roman','FontSize',20);
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
ylim([-1.5 1.1]);

ylabel('Axial stress (Ha/Bohr)','FontSize',22,'FontName','Times New Roman','Interpreter','LaTex');

legend1 = legend('$f_{1}$','$f_{2}$','$f_{3}$','$\sigma$');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

set(figure2,'Units','Inches');
pos = get(figure2,'Position');
set(figure2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'consistency','epsc')
hold off
