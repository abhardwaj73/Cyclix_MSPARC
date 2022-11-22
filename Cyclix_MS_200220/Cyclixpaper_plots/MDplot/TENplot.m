clc
close all
clear all

TEN = dlmread('TEN.txt');

% time
time = [1:1:1300]/1000; % in ps
% shift to 0
time = time - time(1);
TE_NVE = TEN(201:1500);
inset_ind = time < 0.7 & time >= 0.5;
time_inset = time(inset_ind);
% plot NVE total energy of the system
set(0,'defaultTextInterpreter','latex');
set(gca,'XMinorTick','off','YMinorTick','on')
fig1 = figure(1); hold on; box on;
plot(time,TE_NVE,'b','linewidth',1) ;
rectangle('Position',[time_inset(1) -5.8 0.2 0.2])
hold off;
%set(gca,'XTick',[0 2 4 6],'YTick',[-2.75 -2.50 -2.25 -2.00],'TickLength',[0.018 0.025],'FontName','Times  New Roman','FontSize',16);
set(gca,'TickLabelInterpreter','LaTex');
set(gca,'XTick',[0:0.3:time(end)],'YTick',[-6.0:1.0:-3.0],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.1f'),'fontsize',20)
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'),'fontsize',20)
%xlim([0 6.5]);
xlim([0 time(end)]);
%ylim([-2.75 -2.00]);
ylim([-6 -3]);
xlabel('Time (ps)')
%ylabel('system $+$ bath energy (Ha/atom)')
%ylabel('Extended system energy (Ha/atom)')
ylabel('Total energy (Ha/atom)')
set(gca,'FontSize',20);
fig2 = figure(2);
plot(time(inset_ind),TE_NVE(inset_ind),'b','linewidth',1) ;
%set(gca,'XTick',[2.5:0.1:3.0],'YTick',[-2.30705 -2.30700 -2.30695],'TickLength',[0.01 0.015],'FontName','Times New Roman','FontSize',10);
set(gca,'TickLabelInterpreter','LaTex');
set(gca,'XTick',[0.5:0.05:0.7],'YTick',[-5.71178 :0.00004: -5.71166],'TickLength',[0.01 0.015],'FontName','Times New Roman');
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.2f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.5f'))
set(gca,'FontSize',16);
%xlim([0 6.5]);
%xlim([0 time(end)]);
%xlim([time_inset(1) time_inset(end)]);
xlim([0.5 0.7]);
ylim([-5.71178 -5.71166]);

%fig3 = inset(fig1,fig2);
fig3 = inset(fig1,fig2,0.75);
%saveas(gca,'NVTNH_TEX','epsc')

% save picture
h = gcf;
fname = 'NVE_TE';
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
saveas(gcf,fname,'epsc') 
print(h,fname,'-dpdf','-r0')
%print(h,fname,'-dpdf','-r600')
 
 
 
