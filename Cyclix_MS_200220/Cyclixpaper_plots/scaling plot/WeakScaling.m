% plot strong scaling test for (18,9) carbon nanotube
close all; clear all; clc

h = figure;
hold on

nkpt = [9;45;135;270;450;900;1350;2700;4500;11250];%3375;
runtime = [27.893;29.511;32.123;32.291;32.635;33.275;33.835;33.561;34.130;35.760];%34.220;
time_ideal = runtime(1)*ones(size(runtime,1),1);

lg1 = plot(nkpt, runtime, '-o', 'linewidth', 1, 'markersize', 12);
lg1.Color = [0 0 1];

lg2 = plot(nkpt, time_ideal, 'linewidth', 1,'linestyle','--');
lg2.Color = [0 0 0];

set(gca, 'XScale', 'log')
axis([9, 11500, 10, 50]);

set(gca,'XTick',[10 100 1000 10000])
set(gca,'YTick',[10 20 30 40 50])
%ytickformat('%.1f')
set(gca,'YMinorTick','on','TickLength',[0.018 0.025]); % sets the tick mark lengths, major&minor
box on;
hold on;

set(gca,'TickLabelInterpreter','Latex');
set(gca,'fontsize',20,'FontName','Times New Roman');

xlabel('Number of points in Brillouin zone','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Wall time (s)','FontSize',20,'FontName','Times New Roman','Interpreter','LaTex');

legend1 = legend('Cyclix-DFT','Theoretical','Location','Northeast');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'weakscaling','epsc') 
