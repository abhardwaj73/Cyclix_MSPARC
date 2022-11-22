close all
clear
clc

% Convergence of energy, force, stress with eta_points
% System : (12,8) CNT tube, and no external twist
% mesh-size = 0.2 Bohr, eta = 50, C-C bond length = 2.67 Bohr, LDA psp

N_dtp = 7;
n_atm = 2;
kpt = [4 5 6 7 8 9];
kpt_ref = 12;
W_cnt= [-5.6970041963
        -5.7023956807
        -5.7032130931
        -5.7033383134
        -5.7033575669
        -5.7033603342
        %-5.7033607087
        %-5.7033607102
        -5.7033606686
]; % in Ha/atom

F_cnt = [-2.3066292844E-03  -1.3804088796E-02  -1.1567503242E-03
         -7.5900574795E-03   1.1758857135E-02   1.1567503242E-03
         -8.1359769024E-04  -1.4947664161E-02  -1.7889053516E-03
         -6.6719644123E-03   1.3400546953E-02   1.7889053516E-03
         -4.8300180802E-04  -1.5168519900E-02  -1.9198482784E-03
         -6.4560800336E-03   1.3734491988E-02   1.9198482784E-03
         -4.1333185127E-04  -1.5206606140E-02  -1.9486807962E-03
         -6.4070577883E-03   1.3797091925E-02   1.9486807962E-03
         -4.0001066727E-04  -1.5213820363E-02  -1.9521171767E-03
         -6.3978519061E-03   1.3808983652E-02   1.9521171767E-03
         -3.9720681370E-04  -1.5214599237E-02  -1.9531536848E-03
         -6.3956999112E-03   1.3810888723E-02   1.9531536848E-03
         %-3.9694755379E-04  -1.5214850985E-02  -1.9538972467E-03
         %-6.3954718394E-03   1.3811261544E-02   1.9538972467E-03
         % -3.9647905812E-04  -1.5214787734E-02  -1.9540099537E-03
         % -6.3950172479E-03   1.3811353564E-02   1.9540099537E-03
         -3.9663387980E-04  -1.5214950847E-02  -1.9539899342E-03
         -6.3950456916E-03   1.3811318814E-02   1.9539899342E-03
]; % in Ha/Bohr

T_cnt = [7.4034295257E-01
         6.7814534091E-01
         6.6523354297E-01
         6.6321104260E-01 
         6.6281842732E-01
         6.6273018099E-01
         %6.6270940470E-01
         %6.6270379905E-01
         6.6270177710E-01
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
set(gca,'yscale','log','XTick',[4 5 6 7 8 9],'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2],'XMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
%set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',{'$10^{-6}$', '$10^{-5}$','$10^{-4}$', '$10^{-3}$', '$10^{-2}$'});
xlim([3.9 9.1]);
ylim([2.5e-7 1e-2]);

xlabel('Vacuum (Bohr)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Error in energy and force','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');


yyaxis right
plot(kpt,dT_cnt,'b','Marker','d','MarkerSize',12,'LineWidth',1,'LineStyle','--')

set(gca,'TickLabelInterpreter','LaTex');
set(gca,'yscale','log','YTick',[1e-3 1e-2 1e-1 1e0 1e1],'TickLength',[0.018 0.025],'Ycolor',[0 0 0],'FontName','Times New Roman','FontSize',20);
%set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',{'$10^{-3}$','$10^{-2}$','$10^{-1}$', '$10^{0}$','$10^{1}$'});
ylim([3e-3 14]);

ylabel('Error in stress','FontSize',22,'FontName','Times New Roman','Interpreter','LaTex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend1 = legend('Energy (Ha/atom)','Force (Ha/Bohr)','Axial stress (\%)','Location','Northeast');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');


set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'errorVsvac','epsc')
hold off
