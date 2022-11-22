close all
clear
clc

% Convergence of energy, force, stress with mesh size
% System : (12,8) CNT tube, 10 Bohr vacuum, and no external twist
% eta = 50, C-C bond length = 2.67 Bohr, LDA psp

N_dtp = 6;
n_atm = 2;
mesh = [0.30 0.25 0.22 0.175 0.15]';%0.125
mesh_ref = 0.10;
W_cnt= [-5.7009386570
        -5.7025709385
        -5.7032088382
        -5.7033533787
        -5.7033326327
        %-5.7033333899
        -5.7033391452         
]; % in Ha/atom

F_cnt = [-9.6196557199E-04  -2.2070769471E-02   1.7917741842E-02
         -9.7164131792E-03   1.9758673056E-02  -1.7917741842E-02
         -6.5094897510E-04  -1.4910267714E-02  -3.4400271954E-03
         -6.5073375188E-03   1.3431062416E-02   3.4400271954E-03
         -4.1239114593E-04  -1.5302741181E-02  -1.8135486585E-03
         -6.4572079135E-03   1.3869129190E-02   1.8135486585E-03
         -3.5901042025E-04  -1.5248536828E-02  -1.9459410520E-03
         -6.3867285359E-03   1.3849988343E-02   1.9459410520E-03
         -3.5707284697E-04  -1.5279615659E-02  -1.9268981335E-03
         -6.3839645459E-03   1.3886825711E-02   1.9268981335E-03
         % -3.5561976071E-04  -1.5270344640E-02  -1.9371502919E-03
         % -6.3804561165E-03   1.3878217009E-02   1.9371502919E-03
          -3.5689409251E-04  -1.5267546853E-02  -1.9392945300E-03
          -6.3791750207E-03   1.3875583545E-02   1.9392945300E-03
]; % in Ha/Bohr

T_cnt = [6.2229012243E-01
         6.6934809090E-01
         6.6181043091E-01
         6.6158602772E-01
         6.6136042222E-01
         %6.6167294204E-01
         6.6143955479E-01
]; % in Ha/Bohr

% Compute error
%--------------

dW_cnt = abs(W_cnt(1:end-1) - W_cnt(end));
dF_cnt = zeros(N_dtp-1,1);
F_ref = F_cnt(n_atm*(N_dtp-1)+1:end,:);
for i = 1:N_dtp-1
    dF_cnt(i) = 1e-16;
    for j = 1:n_atm
        err = norm(F_cnt(n_atm*(i-1)+j,:) - F_ref(j,:),inf);
        if err > dF_cnt(i)
            dF_cnt(i) = err;
        end
    end
end
dT_cnt = abs((T_cnt(1:end-1) - T_cnt(end))*100/T_cnt(end));


%Linear fit to data
log_dW = log10(dW_cnt);
log_dF = log10(dF_cnt);
log_dT = log10(dT_cnt);
log_mesh = log10(mesh);

P_W = polyfit(log_mesh,log_dW,1);
P_F = polyfit(log_mesh,log_dF,1);
P_T = polyfit(log_mesh,log_dT,1);
W_fit = P_W(1)*log_mesh+P_W(2);
F_fit = P_F(1)*log_mesh+P_F(2);
T_fit = P_T(1)*log_mesh+P_T(2);

antilog_W_fit = 10.^(W_fit);
antilog_F_fit = 10.^(F_fit);
antilog_T_fit = 10.^(T_fit);



% Error plot
%------------

figure1 = figure;
box on

hold on
%S1=scatter(mesh,dW_cnt,144,'r','Marker','o');
S(1)=plot(mesh,dW_cnt,'rd','Marker','o','MarkerSize',12);
S(2)=plot(mesh,dF_cnt,'d','Marker','*','MarkerSize',12);S(2).MarkerEdgeColor = [0 0.6 0];
h(1)=plot(mesh,antilog_W_fit,'r','LineWidth',1,'LineStyle','-');
h(2)=plot(mesh,antilog_F_fit,'color',[0 0.60 0],'LineWidth',1,'LineStyle','-.');
%plot(kpt,dT_cnt,'b','Marker','sq','MarkerSize',10,'LineWidth',1,'LineStyle','-')
set( get( get( h(1), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h(2), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'xscale','log','yscale','log','XTick',[0.15 0.20 0.25 0.30],'XMinorTick','on','xdir','reverse','YTick',[1e-5 1e-4 1e-3 1e-2],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',{'$0.30$','0.25', '$0.20$', '$0.15$'});
set(gca,'YTickLabel',{'$10^{-5}$','$10^{-4}$', '$10^{-3}$', '$10^{-2}$'});
xlim([0.149 0.3015]);
ylim([5e-6 0.022]);

xlabel('Mesh size (Bohr)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Error in energy and force','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

yyaxis right
S(3)=plot(mesh,dT_cnt,'bd','Marker','d','MarkerSize',12);
h(3)=plot(mesh,antilog_T_fit,'b','LineWidth',1,'LineStyle','--');
set( get( get( h(3), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set(gca,'TickLabelInterpreter','LaTex');
set(gca,'yscale','log','YTick',[1e-2 1e-1 1e0 1e1],'TickLength',[0.018 0.025],'Ycolor',[0 0 0],'FontName','Times New Roman','FontSize',20);
%set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',{'$10^{-2}$','$10^{-1}$', '$10^{0}$','$10^{1}$', '$10^{2}$'});
ylim([1e-2 10]);

ylabel('Error in stress','FontSize',22,'FontName','Times New Roman','Interpreter','LaTex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend1 = legend('Energy (Ha/atom)','Force (Ha/Bohr)','Axial stress (\%)','Location','Northeast');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'errorVsmesh','epsc')
hold off
