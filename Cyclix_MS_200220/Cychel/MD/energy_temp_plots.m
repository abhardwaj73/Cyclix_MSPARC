clc
close all
clear all
% NVE results
% TE_NVE = dlmread('TE_NVE');
% TIO_NVE = dlmread('TIO_NVE');
% NVT_NH results
TIO_NVTNH = dlmread('TIO_NVTNH');
TENX_NVTNH = dlmread('TENX_NVTNH') * 27.2114;
% time
time = [1:1:1000];
% plot NVE energy
set(0,'defaultTextInterpreter','latex');

% fig1 = figure(1);
box on

hold on
% plot(time,TE_NVE) ;
% set(gca,'XTick',[0 200 400 600 800 1000],'YTick',[-5.0 -4.0 -3.0 -2.0 -1.0 0.0],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
% set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
% set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
% ylim([-5.0 0.0]);
% xlabel('Time(fs)')
% ylabel('Total energy(Ha/atom)')
% fig2 = figure(2);
% plot(time,TE_NVE) ;
% set(gca,'XTick',[0 500 1000],'YTick',[-3.96194 -3.96190 -3.96186],'TickLength',[0.01 0.015],'FontName','Times New Roman','FontSize',10);
% set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
% set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.5f'))
% ylim([-3.96194 -3.96186]);
% fig3 = inset(fig1,fig2);
% saveas(gca,'NVE_TE','eps')
%title(fig3, 'NVE total energy')
%set(fig3,'Units','Inches');
%pos = get(fig3,'Position');
%set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
%print(fig3,'NVE_TE','-dpdf','-r0')
% hold off
% plot NVE temperature
% fig4 = figure();
% plot(time,TIO_NVE) ;
% set(gca,'XTick',[0 200 400 600 800 1000],'YTick',[100 150 200 250 300 350],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
% set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
% set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.0f'))
% ylim([100 350]);
% xlabel('Time(fs)')
% ylabel('Temperature(K)')
% %title('NVE temperature')
% %saveas(gca,'NVE_TIO','epsc')
% set(fig4,'Units','Inches');
% pos = get(fig4,'Position');
% set(fig4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
% print(fig4,'NVE_TIO','-dpdf','-r0')
% plot NVT energy_extended system
set(0,'defaultTextInterpreter','latex');
fig4 = figure(4);
plot(time,TENX_NVTNH) ;
set(gca,'XTick',[0 250 500 750 1000],'YTick',[-180 -160 -140 -120 -100],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.0f'))
ylim([-180 -100]);
xlabel('Time (fs)')
ylabel('$\mathcal{E}_{system+bath}$ (eV/atom)')
fig5 = figure(5);
plot(time,TENX_NVTNH) ;
set(gca,'XTick',[0 250 500 750 1000],'YTick',[-155.52 -155.51 -155.50],'TickLength',[0.01 0.015],'FontName','Times New Roman','FontSize',10);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.2f'))
ylim([-155.52 -155.50]);
fig6 = inset(fig4,fig5);

%title(fig6, 'NVT extended system energy')
saveas(gca,'NVT_TEX','epsc')

% % plot NVT temperature
set(0,'defaultTextInterpreter','latex');
fig7 = figure(7);
plot(time,TIO_NVTNH) ;
set(gca,'XTick',[0 250 500 750 1000],'YTick',[200 250 300 350 400],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',16);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.0f'))
ylim([200 400]);
xlabel('Time (fs)')
ylabel('Temperature (K)')
fig8 = figure(8);
plot(time,TIO_NVTNH) ;
set(gca,'XTick',[0 250 500 750 1000],'YTick',[299.5 300 300.5],'TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',10);
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.0f'))
set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.1f'))
ylim([299.5 300.5]);
fig9 = inset(fig7,fig8);
%title(fig9,'NVT temperature')
saveas(gca,'NVT_TIO','epsc')

hold off