%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot band structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clc 
% clear

%function BS_plot(fpath,fname,twist,H)
function BS_plot(fpath,fname,Nelectron,H)
format long;
filename = fullfile(fpath,fname);
load(filename);
Ha2eV = 27.211386245988;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Assumed no degeneracy of eigenvalues

%bg_max = 10000;
Nkpoints = size(eign,1);
Nbands = size(eign,2);
%Nbands_occ = ceil(Nelectron/2);
fineness = 1000;
eign_fine = zeros((size(kptsplit_ind,1)-1)*(fineness+1),Nbands);
kpt_fine = zeros((size(kptsplit_ind,1)-1)*(fineness+1),1);
indx_s = 1;

%a = figure;

%a = strcat('figure',ii);
    % a = figure;
    % box on
    % hold on
    % set(gca,'TickLabelInterpreter','LaTex');
    % set(gca,'XTick',[-0.4 -0.2 0.0 0.2 0.4],'XMinorTick','on','YTick',[-10 -8 -6 -4 -2 0 2],'YMinorTick','on','TickLength',[0.018 0.025],'FontName','Times New Roman','FontSize',20);
    % set(gca,'XTickLabel',num2str(get(gca,'XTick')','%.1f'))
    % set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.0f'))
    % xlim([-0.5 0.5]);
    % ylim([-10.5 2]);

for ii = 1:size(kptsplit_ind,1)-1
    [kpt_sort indx_sort] = sort(kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,3));
    eign_1 = eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,:);
    eign_1 = eign_1(indx_sort,:);

    dk = (kpt_sort(end) - kpt_sort(1))/fineness;
    indx_e = indx_s + fineness;
    kpt_fine(indx_s:indx_e,1) = [kpt_sort(1):dk:kpt_sort(end)]';

    for band = 1:Nbands
        eign_fine(indx_s:indx_e,band) = spline(kpt_sort,eign_1(:,band),kpt_fine(indx_s:indx_e,1));
    end
    indx_s = indx_e+1;
    
   

    % y1 = eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,4);
    % y2 = eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,5);

    % [eig_homo, indx_homo] = max(eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,4));
    % [eig_lumo, indx_lumo] = min(eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,5));
    % bg = eig_lumo - eig_homo;
    % if bg < bg_max
    %     bg_max = bg;
    %     kpt_homo = [kpt(kptsplit_ind(ii),2) kpt(kptsplit_ind(ii)+indx_homo-1,3)];
    %     kpt_lumo = [kpt(kptsplit_ind(ii),2) kpt(kptsplit_ind(ii)+indx_lumo-1,3)];
    % end

    % if (ii == 9)
    %     for jj = 4 %1:size(eign,2)
    %         plot(kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,3),eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,jj)*Ha2eV,'b','LineStyle','--','LineWidth',1);
    %         %hasbehavior(p2 ,'legend',false);
    %         hold on;
    %     end

    %     for jj = 5  %1:size(eign,2)
    %         CB = plot(kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,3),eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,jj)*Ha2eV,'b','LineStyle','--','LineWidth',1);
    %         set( get( get( CB, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    %         hold on;
    %     end
    % end

    % if (ii == 4)
    %     for jj = 4 %1:size(eign,2)
    %         plot(kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,3),eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,jj)*Ha2eV,'r','LineStyle','-','LineWidth',1);
    %         %hasbehavior(p2 ,'legend',false);
    %         hold on;
    %     end

    %     for jj = 5  %1:size(eign,2)
    %         CB = plot(kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,3),eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,jj)*Ha2eV,'r','LineStyle','-','LineWidth',1);
    %         set( get( get( CB, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    %         hold on;
    %     end
    % end

    
    % legend(axes1,'show');
    % set(legend,'Position',[0.701755847995657 0.811523808621225 0.209583544871275 0.110285717510042],'Interpreter','latex');
    % set(legend,'Position',[0.511874832391208 0.859333333236831 0.400059861794459 0.0622857158978779],'Interpreter','latex','Orientation','horizontal');
    % xlabel('$\frac{\eta H} {2\pi}$','Interpreter','latex','FontSize',20,'FontName','Times New Roman');
    % ylabel('Energy (eV)','Interpreter','latex','FontSize',20,'FontName','Times New Roman');
    % set(a,'Units','Inches');
    % pos = get(a,'Position');
    % %disp([1.02*pos(3), 1.02*pos(4)]);aa
    % set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
    % filename = strcat('C_',num2str(ii));
    % file_loc = fullfile(fpath,filename);
    % saveas(gcf,file_loc,'epsc')
    %hold off
    
end

% xlabel('$\frac{\eta H} {2\pi}$','Interpreter','latex','FontSize',20,'FontName','Times New Roman');
% ylabel('Energy (eV)','Interpreter','latex','FontSize',20,'FontName','Times New Roman');


% set(a,'Units','Inches');
% %legend()
% pos = get(a,'Position');
%     %disp([1.02*pos(3), 1.02*pos(4)]);aa
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4,2])
% filename = strcat('C_',num2str(ii));
% file_loc = fullfile(fpath,filename);
% saveas(gcf,file_loc,'epsc')

% hold off

% eign_split1 = [eign(:,2:end) eign(:,1)];
% eign_split2 = [eign(:,end) eign(:,1:end-1)];
% %eign_cur = eign(kptsplit_ind(ii):kptsplit_ind(ii+1)-1,:); 

% band_homo = 1./(1+ exp((eign - E_f)/(1e-6))) > 1e-3 & 1./(1+ exp((eign_split1 - E_f)/(1e-6))) < 1e-3;
% band_lumo = 1./(1+ exp((eign - E_f)/(1e-6))) < 1e-3 & 1./(1+ exp((eign_split2 - E_f)/(1e-6))) > 1e-3;

% x = kpt(kptsplit_ind(ii):kptsplit_ind(ii+1)-2,3);
% y1 = eign_cur(band_homo);
% y2 = eign_cur(band_lumo);

% y1 = y1(1:end-1);
% y2 = y2(1:end-1);

% dx = (kpt(kptsplit_ind(ii+1)-2,3) - kpt(kptsplit_ind(ii),3))/1000;
% xx = [kpt(kptsplit_ind(ii),3):dx:kpt(kptsplit_ind(ii+1)-2,3)];
% yy1 = spline(x,y1,xx);
% yy2 = spline(x,y2,xx);

% [eig_homo, indx_homo] = max(yy1);
% [eig_lumo, indx_lumo] = min(yy2);
% bg = eig_lumo - eig_homo;
% if bg < bg_max
%     bg_max = bg;
%     kpt_homo = [kpt(kptsplit_ind(ii),2) xx(indx_homo)];
%     kpt_lumo = [kpt(kptsplit_ind(ii),2) xx(indx_lumo)];
% end




eign_fine_all = reshape(eign_fine',[],1);
[eign_fine_sort,indx_fine_sort] = sort(eign_fine_all);

Nkpoints_fine = size(eign_fine,1);
Nbands_fine_occ = ceil(Nelectron*Nkpoints_fine/2);


bg = eign_fine_sort(Nbands_fine_occ+1) - eign_fine_sort(Nbands_fine_occ);
indx_homo = ceil(indx_fine_sort(Nbands_fine_occ)/Nbands);
indx_lumo = ceil(indx_fine_sort(Nbands_fine_occ+1)/Nbands);
kpt_homo = [kpt(kptsplit_ind(ceil(indx_homo/(fineness+1))),2) kpt_fine(indx_homo)];
kpt_lumo = [kpt(kptsplit_ind(ceil(indx_lumo/(fineness+1))),2) kpt_fine(indx_lumo)];

band_homo = rem(indx_fine_sort(Nbands_fine_occ),Nbands);
band_lumo = rem(indx_fine_sort(Nbands_fine_occ+1),Nbands); % assumed number of bands are greater than Nbands_occ+1

fprintf("With spline: Band gap = %.10f\n kpt_homo = [%f %f]\n kpt_lumo = [%f %f]\n band_homo %d band_lumo %d\n",bg,kpt_homo(1),kpt_homo(2),kpt_lumo(1),kpt_lumo(2),band_homo,band_lumo);


% Without spline

eign_all = reshape(eign',[],1);
[eign_sort,indx_sort] = sort(eign_all);

Nbands_occ = ceil(Nelectron*Nkpoints/2);
bg = eign_sort(Nbands_occ+1) - eign_sort(Nbands_occ);
indx_homo = ceil(indx_sort(Nbands_occ)/Nbands);
indx_lumo = ceil(indx_sort(Nbands_occ+1)/Nbands);

kpt_homo = [kpt(indx_homo,2) kpt(indx_homo,3)];
kpt_lumo = [kpt(indx_lumo,2) kpt(indx_lumo,3)];

band_homo = rem(indx_sort(Nbands_occ),Nbands);
band_lumo = rem(indx_sort(Nbands_occ+1),Nbands);

%eign_sort(Nbands_occ-5:Nbands_occ+5);

% [eig_homo, indx_homo] = max(eign(:,4));
% [eig_lumo, indx_lumo] = min(eign(:,5));
% bg = eig_lumo - eig_homo;
% kpt_homo = [kpt(indx_homo,2) kpt(indx_homo,3)];
% kpt_lumo = [kpt(indx_lumo,2) kpt(indx_lumo,3)];

fprintf("Without spline: Band gap = %.10f\n kpt_homo = [%f %f]\n kpt_lumo = [%f %f]\nband_homo %d band_lumo %d",bg,kpt_homo(1),kpt_homo(2),kpt_lumo(1),kpt_lumo(2),band_homo,band_lumo);

% Calculate fermi velocity and effective mass of electrons and holes
%-------------------------------------------------------------------

fac = (2*pi/H);
kpt(:,3) = kpt(:,3)*fac;

% mass of electron
h_bar = 6.62607015 * (1e-34)/(2*pi);
J2eV = 6.241509074 * (1e18);
Bohr2m = 5.29177210903 * 1e-11;
m_e_0 = 9.1093837015 * 1e-31;

m0 = (h_bar^2 * J2eV / (Bohr2m^2))/m_e_0;

% Hole's effective mass
for ii = 1:size(kptsplit_ind,1)
    if kpt(kptsplit_ind(ii),2) == kpt_homo(1)
        indx_nu_homo = ii;
        break;
    end
end 

range_eta_homo = kptsplit_ind(indx_nu_homo):kptsplit_ind(indx_nu_homo+1)-1;

homo_spline = spline(kpt(range_eta_homo,3),eign(range_eta_homo,band_homo));
%v_h = ppval(fnder(homo_spline,1),kpt_homo(2)*fac);
m_h = 1/ppval(fnder(homo_spline,2),kpt_homo(2)*fac);

V_h = ppval(fnder(homo_spline,1),kpt_fine(1:1000)*fac);
M_h = 1./ppval(fnder(homo_spline,2),kpt_fine(1:1000)*fac);

M_h(abs(V_h) < 1e-2)

% Electron's effective mass
for ii = 1:size(kptsplit_ind,1)
    if kpt(kptsplit_ind(ii),2) == kpt_lumo(1)
        indx_nu_lumo = ii;
        break;
    end
end

range_eta_lumo = kptsplit_ind(indx_nu_lumo):kptsplit_ind(indx_nu_lumo+1)-1;

lumo_spline = spline(kpt(range_eta_lumo,3),eign(range_eta_lumo,band_lumo));
%v_e = ppval(fnder(lumo_spline,1),kpt_lumo(2)*fac);
m_e = 1/ppval(fnder(lumo_spline,2),kpt_lumo(2)*fac);


V_e = ppval(fnder(lumo_spline,1),kpt_fine(1:1000)*fac);
M_e = 1./ppval(fnder(lumo_spline,2),kpt_fine(1:1000)*fac);


M_e(abs(V_e) < 1e-2)


fprintf("\n\nEffective mass: Hole = %f, Electron = %f\n",m_h,m_e);

% True Fermi velocity (valid only for metals)
%--------------------------------------------

if(bg < 1e-3)

    dp = 10;
    l = -1;

    v0 = 2.18769126364 * 1e6;

    v_h = max(abs((eign(indx_homo:indx_homo+dp,band_homo) - eign(indx_homo+l:indx_homo+l+dp,band_homo))./(kpt(indx_homo:indx_homo+dp,3)-kpt(indx_homo+l:indx_homo+l+dp,3))))*v0
    v_e = max(abs((eign(indx_lumo:indx_lumo+dp,band_lumo) - eign(indx_lumo+l:indx_lumo+l+dp,band_lumo))./(kpt(indx_lumo:indx_lumo+dp,3)-kpt(indx_lumo+l:indx_lumo+l+dp,3))))*v0

    fprintf("\n\nFermi velocity: Hole = %f, Electron = %f\n",v_h,v_e);
end