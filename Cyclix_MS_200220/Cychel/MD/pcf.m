function Pcf = pcf(equilStep, stepPcf, limitPcf, sizeHistPcf, rangePcf, cell, pos, n_types, typ_natm)
% Assumption: atoms of different types are segregated and provided as an array one after the other

% Pair Correlation Function =>
% Input variables:
% equilStep = time step at which equilibration is achieved (will collect
% data after this step)
% stepPcf = intervals at which data is to be used to plot pcf
% limitPcf = total number of time steps to be used for pcf
% histPcf = count the number of atoms in a spherical shell of thickness
% deltaR
% sizeHistPcf = total number of spherical shells to be considered
% rangePcf = maximum radial distance to be  considered 
% cell = vector comprising of length of unit cell
% pos = an array of size [tn_atoms*n_ts,3] containing atomic positions for all time steps
% n_types = number of atom types
% typ_natm = an array containing number of atoms of each type

% Output variables:
% Pcf: pair correlation function for different radial distances

tn_atoms = sum(typ_natm);
n_ts = size(pos,1)/tn_atoms;
last = stepPcf * limitPcf + equilStep;
if last > n_ts
    error('Error: "%d time steps for calulating PCF @ steps of %d are unavailable!" \n',limitPcf,stepPcf);
end

Traj = zeros(tn_atoms,3,n_ts);
for i = 1:n_ts
    Traj(:,:,i) = pos((i-1)*tn_atoms+1:i*tn_atoms,:);
end
atm_filter = Traj(:,:,(equilStep+1:stepPcf:last));
deltaR = rangePcf/sizeHistPcf;
npcf = n_types*n_types;
histPcf = zeros(sizeHistPcf,npcf);
Pcf = zeros(sizeHistPcf,npcf);
seq = [0.5:sizeHistPcf-0.5]' ;
volume = prod(cell) ;

count = 0;
for typ = 1:n_types
  for j1 = count+1:count+typ_natm(typ)
    for j2 = j1+1:count+typ_natm(typ)
      dr = abs(atm_filter(j1,:,:) - atm_filter(j2,:,:));
      dr(1,1,dr(1,1,:)>0.5*cell(1)) = abs(dr(1,1,dr(1,1,:)>0.5*cell(1)) - cell(1)) ;
      dr(1,2,dr(1,2,:)>0.5*cell(2)) = abs(dr(1,2,dr(1,2,:)>0.5*cell(2)) - cell(2)) ;
      dr(1,3,dr(1,3,:)>0.5*cell(3)) = abs(dr(1,3,dr(1,3,:)>0.5*cell(3)) - cell(3)) ;
      rr = sqrt((dr(1,1,:)).^2 + (dr(1,2,:)).^2 + (dr(1,3,:)).^2) ;
      rr = reshape(rr,[],1,1) ;        
      n = ceil((rr(rr < rangePcf))/deltaR) ; 
      TF = isempty(n);       
      if TF == 0
        [element,~,~]  = unique(n) ;
        %times = histcounts(cc,'BinMethod','integers');
        times = histc(n,unique(n));
        histPcf(element,typ) = histPcf(element,typ) + times ;
      end
    end
  end
  count = count + typ_natm(typ);
  normfac = volume/(2*pi*(deltaR^3)*typ_natm(typ)*typ_natm(typ)*limitPcf) ;
  Pcf(:,typ) = (histPcf(:,typ) * normfac)./(seq.^2) ;
end
cc = n_types+1;
count = 0;
for typ1 = 1:n_types
  count2 = 0;
  for typ2 = 1:n_types
    if typ1 == typ2
    else
      for j1 = count+1:count+typ_natm(typ1)
        for j2 = count2 + 1: count2 + typ_natm(typ2) 
          dr = abs(atm_filter(j1,:,:) - atm_filter(j2,:,:));
          dr(1,1,dr(1,1,:)>0.5*cell(1)) = abs(dr(1,1,dr(1,1,:)>0.5*cell(1)) - cell(1)) ;
          dr(1,2,dr(1,2,:)>0.5*cell(2)) = abs(dr(1,2,dr(1,2,:)>0.5*cell(2)) - cell(2)) ;
          dr(1,3,dr(1,3,:)>0.5*cell(3)) = abs(dr(1,3,dr(1,3,:)>0.5*cell(3)) - cell(3)) ;
          rr = sqrt((dr(1,1,:)).^2 + (dr(1,2,:)).^2 + (dr(1,3,:)).^2) ;
          rr = reshape(rr,[],1,1) ;        
          n = ceil((rr(rr < rangePcf))/deltaR) ; 
          TF = isempty(n);       
          if TF == 0
            [element,~,~]  = unique(n) ;
            %times = histcounts(cc,'BinMethod','integers');
            times = histc(n,unique(n));
            histPcf(element,cc) = histPcf(element,cc) + times ;
          end
        end
      end
      normfac = volume/(2*pi*(deltaR^3)*(typ_natm(typ1))^2*limitPcf) ;
      Pcf(:,cc) = (histPcf(:,cc) * normfac)./(seq.^2) ;
      cc = cc + 1;
    end
    count2 = count2 + typ_natm(typ2);
  end
  count = count + typ_natm(typ1);
end
radial_grid = 0:deltaR:rangePcf-deltaR ;
% Only for two atoms
if(Pcf(:,3) ~= Pcf(:,4))
  fprintf('Warning: correlation 1-2 is not same as 2-1\n');
end

fileID = fopen('Out_pcf.txt','w');
for ii = 1:n_types
  fprintf(fileID,':R: \t :pcf(%d -- %d):\n',ii,ii);
  for kk = 1:sizeHistPcf
    fprintf(fileID,'%.15f %.15f\n',radial_grid(kk),Pcf(kk,ii));
  end
end

cc = n_types + 1;
for typ1 = 1:n_types
  for typ2 = typ1+1:n_types
    fprintf(fileID,':R: \t :pcf(%d -- %d):\n',typ1,typ2);
    for kk = 1:sizeHistPcf
      fprintf(fileID,'%.15f %.15f\n',radial_grid(kk),Pcf(kk,cc));
    end
    cc = cc + 1;
  end
  cc = cc + typ1;
end

fclose(fileID);
% Plots:

% for ii = 1:n_types
%   figure()
%   plot(radial_grid,Pcf(:,ii)) ;
%   xlabel('Radial distance from the reference atom (in Bohr)')
%   ylabel('g(r)')
%   title( ['PCF for ',num2str(ii),'--',num2str(ii)])
%   saveas(gcf,sprintf('FIG%d',ii),'epsc')
% end
% count = n_types + 1;
% for ii = 1:n_types
%   for jj = 1:n_types
%     if ii == jj
%       continue
%     else
%       figure()
%       plot(radial_grid,Pcf(:,count)) ;
%       xlabel('Radial distance from the reference atom (in Bohr)')
%       ylabel('g(r)')
%       title( ['PCF for ',num2str(ii),'--',num2str(jj)])
%       saveas(gcf,sprintf('FIG%d',count),'epsc')
%       count = count + 1;
%     end
%   end
% end