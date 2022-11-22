function S = ReadPseudoPot(lloc,element)
% @brief ReadPseudoPot(lloc,element,dx,dy,dz) makes the Pseudopotential files
% @param lloc		Index of angular momentum chosen to be the local components of 
%					the pseudopotential
% @param element	Atom type name
% @param dx			Grid size in x-direction
% @param dy			Grid size in y-direction
% @param dz			Grid size in z-direction
%
% @authors	Qimen Xu <qimenxu@gatech.edu>
%			Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% @2016 Georgia Institute of Technology.
%

%filename = strcat('./Pseudopotentials/psd_',element,'.pot');
filename = sprintf('%s/%s',S.inputfile_path, psdfname);
[fid,~] = fopen(filename,'r');  %file opened
assert(fid~=-1,'Error: Check the element name and the corresponding .pot file') ;
fscanf(fid,'%s',1);  % Name
str = fgetl(fid);
while(strcmp(str,' Radial grid follows') == false)
    str = fgetl(fid);
end

% radial grid
r = fscanf(fid,'%g'); 

% read pseudopotential (l = 0,1,2,3)
str = fgetl(fid);
while(strcmp(str,'Pseudopotential follows (l on next line)') == true)
    l = fscanf(fid,'%g',1);
	if ((l - round(l)) ~= 0)
		error('Azimuthal quantum number <l> must be an interger, check .pot file');
	elseif (l<0)||(l>3)
		error('Only <l> between 0 and 3 supported');
	end	
	V(:,l+1) = fscanf(fid,'%g');
    str = fgetl(fid);
end
lmax = l;
U = zeros(length(r),lmax+1);
rc_v = zeros(lmax+1,1);

while(strcmp(str,'Core charge follows') == true)
    CoreCharge = fscanf(fid,'%g'); % CoreCharge
    str = fgetl(fid);
end
while(strcmp(str,'Valence charge follows') == true)
    uu = fscanf(fid,'%g');
    str = fgetl(fid);
end
while(strcmp(str,'Pseudo-wave-function follows (l, zelect, rc)') == true)
    l = fscanf(fid,'%g',1);
	fscanf(fid,'%g',1); % Z
	rc_v(l+1) = fscanf(fid,'%g',1);
	U(:,l+1) = fscanf(fid,'%g');	
    str = fgetl(fid);
end
fclose(fid);

% Conversion from Rydberg to Hartree, Division by r and r^2
V = 0.5 * V ./ (r*ones(1,lmax+1));
U = U ./ (r*ones(1,lmax+1));

uu = uu ./ (4*pi*r.^2);

r_guess = r;
rho_isolated_guess = uu;
% Including the point r=0 in the pseudopotential data
r=[0;r];
V = [V(1,:);V];
U = [U(1,:);U];
r_guess = [0;r_guess];
rho_isolated_guess = [rho_isolated_guess(1);rho_isolated_guess];

% setting local component
Vloc = V(:,lloc+1);

UdV = zeros(length(r),lmax+1);
for k = 0:lmax
	UdV(:,k+1) = U(:,k+1).*(V(:,k+1)-Vloc);
end


% Calculate denominator for non-local
rc = max(rc_v); % cut-off radius
dr = 1e-6;
r_temp = 0:dr:rc+1; % integrate over 0 to rc + 1 bhor
w = dr * ones(size(r_temp)); w(1)=0.5;w(end)=0.5;
% gamma_Jl = ones(lmax+1,1);
Denom = ones(lmax+1,1);
for k = 0:lmax
	UdV_temp = interp1(r, UdV(:,k+1), r_temp, 'spline');
	UJl_temp = interp1(r, U(:,k+1), r_temp, 'spline');
	if k ~= lloc
		Denom(k+1,1) = sum(UdV_temp.*UJl_temp.*w.*r_temp.*r_temp);
	end
end
gamma_Jl = 1./Denom;
% gamma_Jl = 0*gamma_Jl;

% Saving the pseudopotential data to matlab .mat files
r_grid_vloc = r;
r_grid_rho = r;
% PATH = strcat('./Pseudopotentials/Pseudopotential_', element);
% save(PATH,'Vloc','r_grid_vloc','rc','gamma_Jl','UdV','lmax') ;
% PATH = strcat('./Pseudopotentials/Isolated_atom_rho_',element);
% save(PATH,'r_grid_rho','rho_isolated_guess');

S.Atm(ityp).Z = Z;

S.Atm(ityp).Vloc = Vloc;
S.Atm(ityp).r_grid_vloc = r_grid_vloc;
S.Atm(ityp).rc = rc; 
S.Atm(ityp).gamma_Jl = gamma_Jl;
S.Atm(ityp).UdV = UdV;
S.Atm(ityp).lmax = lmax;
S.Atm(ityp).lloc = lloc;
S.Atm(ityp).r_grid_rho = r_grid_rho;
S.Atm(ityp).rho_isolated_guess = rho_isolated_guess;
