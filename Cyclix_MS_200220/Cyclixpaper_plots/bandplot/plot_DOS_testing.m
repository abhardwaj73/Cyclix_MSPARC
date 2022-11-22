close all;

% set up lambda here
% TODO: 
filename = fullfile('./18_9/0','SPARCEigenvalues_C.mat');
load(filename);
%lambda = eign;


% Fineness of k-grid
fineness = 1000;
Nkpoints = size(eign,1);
Nbands = size(eign,2);
eign_fine = zeros((1*(fineness+1)),2);%Nbands);%size(kptsplit_ind,1)-1)
kpt_fine = zeros((1*(fineness+1)),1);%size(kptsplit_ind,1)-1)
indx_s = 1;

for ii = 4 %1:size(kptsplit_ind,1)-1
	[kpt_sort indx_sort] = sort(kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,3));
	eign_1 = eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,4:5);
	eign_1 = eign_1(indx_sort,:);

	dk = (kpt_sort(end) - kpt_sort(1))/fineness;
	indx_e = indx_s + fineness;
	kpt_fine(indx_s:indx_e,1) = [kpt_sort(1):dk:kpt_sort(end)]';

	for band = 1:2
	    eign_fine(indx_s:indx_e,band) = spline(kpt_sort,eign_1(:,band),kpt_fine(indx_s:indx_e,1));
	end
	indx_s = indx_e+1;
end


lambda = eign_fine;



% plot Density Of States (DOS)
eV2Ha = 1 / 27.21138397;

% mu = 0.03 * eV2Ha;
% N = 1 * size(lambda,1) * size(lambda,2);
% [DOS,E] = eig2DOS(lambda, N, mu);

% h = figure;
% plot(E, DOS);
% xlim([min(min(lambda)), max(max(lambda))]);

% xlabel('$E$ (Ha)', 'interpreter', 'latex');
% %ylabel('Density of states/Ha', 'interpreter', 'latex');
% ylabel('Density of states', 'interpreter', 'latex');

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
% saveas(gcf,'DOS_Ha','epsc') 
% print(h,'DOS_Ha','-dpdf','-r0')


% Plot DOS in eV
% convert lambda to eV
lambda_eV = lambda / eV2Ha;
mu = 0.027; % in eV
N = size(lambda_eV,1) * size(lambda_eV,2);
[DOS_eV,E_eV] = eig2DOS(lambda_eV, N, mu);

h_eV = figure;
box on
hold on
plot(DOS_eV,E_eV,'r','LineStyle','-','LineWidth',1);


filename = fullfile('./18_9','SPARCEigenvalues_C.mat');
load(filename);
%lambda = eign;


%Fineness of k-grid
fineness = 1000;
Nkpoints = size(eign,1);
Nbands = size(eign,2);
eign_fine = zeros((1*(fineness+1)),2);%Nbands);%size(kptsplit_ind,1)-1)
kpt_fine = zeros((1*(fineness+1)),1);%size(kptsplit_ind,1)-1)
indx_s = 1;

for ii = 4
	[kpt_sort indx_sort] = sort(kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,3));
	eign_1 = eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,4:5);
	eign_1 = eign_1(indx_sort,:);

	dk = (kpt_sort(end) - kpt_sort(1))/fineness;
	indx_e = indx_s + fineness;
	kpt_fine(indx_s:indx_e,1) = [kpt_sort(1):dk:kpt_sort(end)]';

	for band = 1:2 % 1:Nbands
	    eign_fine(indx_s:indx_e,band) = spline(kpt_sort,eign_1(:,band),kpt_fine(indx_s:indx_e,1));
	end
	indx_s = indx_e+1;
end


lambda = eign_fine;



lambda_eV = lambda / eV2Ha;
mu = 0.027; % in eV
N = size(lambda_eV,1) * size(lambda_eV,2);
[DOS_eV,E_eV] = eig2DOS(lambda_eV, N, mu);



plot(DOS_eV,E_eV,'b','LineStyle','--','LineWidth',1);




















set(gca,'TickLabelInterpreter','LaTex');
set(gca,'YTick',[-10 -8 -6 -4 -2 0 2],'YMinorTick','on','XTick',[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7],'XMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.1f'))

ylim([-10.5,2]);
xlim([0.0,0.7]);

%ylabel('Energy (eV)','Interpreter','latex','FontSize',20,'FontName','Times New Roman');
xlabel('Density of states (eV$^{-1}$)','Interpreter','latex','FontSize',20,'FontName','Times New Roman');

legend1 = legend('$\gamma = 0.0\%$','$\gamma = 2.8\%$','Location','best');
set(legend1,'fontsize',20,'FontName','Times New Roman','Interpreter','LaTex');

set(h_eV,'Units','Inches');
pos = get(h_eV,'Position');
set(h_eV,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
saveas(gcf,'DOS_eV','epsc')
