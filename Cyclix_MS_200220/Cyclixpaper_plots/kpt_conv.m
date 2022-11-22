close all
clear
clc

% Convergence of energy, force, stress with eta_points
% System : (18,9) CNT tube, 10 Bohr vacuum, and no external twist
% mesh-size = 0.2 Bohr, C-C bond length = 2.67 Bohr, LDA psp

N_dtp = 11;
n_atm = 2;
kpt = [1 20 40 60 100 120 140 160 180 200];%80 
kpt_ref = 300;
W_cnt= [
        -5.6700264121
        -5.7033890786
        -5.7033508652
        -5.7033575300
        %-5.7033615031
        -5.7033620316
        -5.7033616867
        -5.7033614546
        -5.7033614423
        -5.7033614763
        -5.7033615022
        -5.7033614895        
]; % in Ha/atom



F_cnt = [
         -4.4859123566E-03  -3.1879412221E-02  -1.0904113814E-01
         -1.6754447139E-02   2.7490462339E-02   1.0904113814E-01
         -3.5850411250E-04  -1.6708734683E-02  -2.5166530924E-03
         -6.9524594428E-03   1.5197935996E-02   2.5166530924E-03
         -4.1680576147E-04  -1.5199653459E-02  -2.4996217543E-03
         -6.4076753271E-03   1.3789315731E-02   2.4996217543E-03
         -4.3638426864E-04  -1.4574650087E-02  -2.3579755768E-03
         -6.1779055135E-03   1.3207791353E-02   2.3579755768E-03
         % -4.2681958619E-04  -1.4650021067E-02  -2.1560591346E-03
         % -6.1989234233E-03   1.3280635642E-02   2.1560591346E-03
         -4.2175849701E-04  -1.4752021978E-02  -2.1213174611E-03
         -6.2347407283E-03   1.3376417902E-02   2.1213174611E-03
         -4.2161727580E-04  -1.4782762490E-02  -2.1474713219E-03
         -6.2468188470E-03   1.3404693819E-02   2.1474713219E-03
         -4.2332963054E-04  -1.4767866246E-02  -2.1732885934E-03
         -6.2424134807E-03   1.3390339120E-02   2.1732885934E-03
         -4.2392177302E-04  -1.4751209481E-02  -2.1766819195E-03
         -6.2364101882E-03   1.3374844838E-02   2.1766819195E-03
         -4.2367876885E-04  -1.4748838355E-02  -2.1697434427E-03
         -6.2352687843E-03   1.3372738637E-02   2.1697434427E-03
         -4.2328856401E-04  -1.4753075444E-02  -2.1658894013E-03
         -6.2365533629E-03   1.3376779009E-02   2.1658894013E-03
         -4.2338857317E-04  -1.4754226713E-02  -2.1676783039E-03
         -6.2371435465E-03   1.3377800071E-02   2.1676783039E-03
]; % in Ha/Bohr

T_cnt = [ 
        -1.3859052176E-01 
         6.7119696261E-01
         6.5842529964E-01
         6.5444439856E-01
         %6.5647386230E-01
         6.5753558597E-01
         6.5759896573E-01
         6.5729902909E-01
         6.5713994838E-01
         6.5717078235E-01
         6.5723160722E-01
         6.5722835765E-01
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

% Error plot
%------------

figure1 = figure;
box on

hold on
plot(kpt,dW_cnt,'r','Marker','o','MarkerSize',12,'LineWidth',1,'LineStyle','-')
plot(kpt,dF_cnt,'color',[0 0.60 0],'Marker','*','MarkerSize',12,'LineWidth',1,'LineStyle','-.')
%plot(kpt,dT_cnt,'b','Marker','sq','MarkerSize',10,'LineWidth',1,'LineStyle','-')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'yscale','log','XTick',[0 40 80 120 160 200],'YTick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'XMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
%set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',{'$10^{-8}$','$10^{-7}$','$10^{-6}$', '$10^{-5}$','$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
xlim([0 202]);
ylim([1e-8 0.2]);

xlabel('Number of $\eta$-points','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Error in energy and force','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');


yyaxis right
plot(kpt,dT_cnt,'b','Marker','d','MarkerSize',12,'LineWidth',1,'LineStyle','--')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'yscale','log','YTick',[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'TickLength',[0.018 0.025],'Ycolor',[0 0 0],'FontName','Times New Roman','FontSize',20);
%set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$', '$10^{0}$','$10^{1}$', '$10^{2}$'});
ylim([1e-4 150]);

ylabel('Error in stress','FontSize',22,'FontName','Times New Roman','Interpreter','LaTex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend1 = legend('Energy (Ha/atom)','Force (Ha/Bohr)','Axial stress (\%)','Location','Northeast');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');


set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'errorVseta','epsc')
hold off
