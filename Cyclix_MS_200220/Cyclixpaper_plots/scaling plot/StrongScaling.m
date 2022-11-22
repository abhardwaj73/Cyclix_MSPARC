% plot strong scaling test for (18,9) carbon nanotube
close all; clear all; clc

h = figure;
hold on

np = [1;5;10;15;30;45;54;90;135;270;540;1080];%18;27;
runtime = [2161.161;484.918;260.303;172.014;88.865;61.389;51.207;32.435;21.730;12.636;6.927;4.035]; %149.865;98.520;
time_ideal = runtime(1)*np(1)./np;

lg1 = loglog(np, runtime, '-o', 'linewidth', 1, 'markersize', 12);
lg1.Color = [0 0 1];

lg2 = loglog(np, time_ideal, 'linewidth', 1,'linestyle','--');
lg2.Color = [0 0 0];

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis([1, 1100, 1, 2200]);

set(gca,'XTick',[1 10 100 1000])
%ytickformat('%.1f')
set(gca,'TickLength',[0.018 0.025]); % sets the tick mark lengths, major&minor
box on;
hold on;

set(gca,'TickLabelInterpreter','Latex');
set(gca,'fontsize',20,'FontName','Times New Roman');


xlabel('Number of CPU cores','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Wall time (s)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

legend1 = legend('Cyclix-DFT','Theoretical','Location','Northeast');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'strongscaling','epsc') 
