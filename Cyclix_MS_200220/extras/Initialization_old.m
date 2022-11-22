function S = Initialization(filename)
% @brief	Initialization(filename) reads data from input file 'filename'
%			and creates a struct 'S' to store all the data read and calculated.
%
% @param filename	The data filename
%
% @authors	Qimen Xu <qimenxu@gatech.edu>
%			Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @2016 Georgia Institute of Technology.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Read Input file                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Reading Input file...\n');
fid1=fopen(filename,'r');
textscan(fid1,'%s',1,'delimiter','\n');
C_BC = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_meth = textscan(fid1,'%d',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_vec1 = textscan(fid1,'%f %f %f',1,'delimiter','\n');
C_vec2 = textscan(fid1,'%f %f %f',1,'delimiter','\n');
C_vec3 = textscan(fid1,'%f %f %f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
sides = textscan(fid1,'%f %f %f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_int = textscan(fid1,'%f %f %f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_kpt = textscan(fid1,'%f %f %f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_trevsym = textscan(fid1,'%d',1,'delimiter','\n');
% textscan(fid1,'%s',2,'delimiter','\n');
% C_bs = textscan(fid1,'%d',1,'delimiter','\n');
% textscan(fid1,'%s',2,'delimiter','\n');
% C_lattice = textscan(fid1,'%s',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_fdn = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_npl = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_conditioner = textscan(fid1,'%d',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_thet = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_maxrlx = textscan(fid1,'%d',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_scf = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_nev = textscan(fid1,'%d',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_poist = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_pseudchgt = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_rcref = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_typ = textscan(fid1,'%f',1,'delimiter','\n');
textscan(fid1,'%s',2,'delimiter','\n');
C_dbg = textscan(fid1,'%f',1,'delimiter','\n');
dbg_switch = C_dbg{1,1}; % Debug switch
textscan(fid1,'%s',2,'delimiter','\n');
C_nspin = textscan(fid1,'%f',1,'delimiter','\n');
nspin = C_nspin{1,1};
textscan(fid1,'%s',2,'delimiter','\n');
C_xc = textscan(fid1,'%d',1,'delimiter','\n');
xc = C_xc{1,1};

n_typ = C_typ{1,1}; % no. of atom types
n_atm=0; Nelectron = 0; Atoms = [];
% Atom types info stored in structure Atm with fields typ, Z, coords
for ii=1:n_typ
    textscan(fid1,'%s',2,'delimiter','\n');
    C_sym = textscan(fid1,'%s',1,'delimiter','\n');
    C_psptyp = textscan(fid1,'%d',1,'delimiter','\n');
    C_chg = textscan(fid1,'%f',1,'delimiter','\n');
    C_mass = textscan(fid1,'%f',1,'delimiter','\n');
    C_natm = textscan(fid1,'%f',1,'delimiter','\n');
    textscan(fid1,'%s',1,'delimiter','\n');
    C_lloc = textscan(fid1,'%d',1,'delimiter','\n');
    % local component of the pseudopotential - 0,1,2 or 3, for s,p,d & f respectively
    lloc = C_lloc{1,1};
    textscan(fid1,'%s',1,'delimiter','\n');
    C_atom = textscan(fid1,'%f %f %f',C_natm{1,1},'delimiter','\n');
    textscan(fid1,'%s',1,'delimiter','\n');
    C_mag = textscan(fid1,'%f',C_natm{1,1},'delimiter','\n');
    if (feof(fid1)~=1 && ii==n_typ) || (feof(fid1)==1 && ii~=n_typ)
        error('Error: Inconsistency between no. of types of atoms and the ones actually specified!')
    end
    typ = C_sym{1,1}; % element type
    psptyp = C_psptyp{1,1};
    Z = C_chg{1,1}; % Charge on nucleus
    mass = C_mass{1,1}; % mass of atom
    n_atm_typ = C_natm{1,1}; % no. of atoms of each type
    n_atm = n_atm + C_natm{1,1};
    Nelectron = Nelectron + Z*n_atm_typ;
    coords = [C_atom{:,1},C_atom{:,2},C_atom{:,3}];
    mag = C_mag{:};
    if(sum(abs(mag(:)) > Z))
        error("Magenetization specified for atom typ %d exceeds number of valence electrons!\n",ii);
    end
    Atoms = [Atoms;coords]; % full list of atom coords of all types
    Atm(ii)=struct('typ',typ,'psptyp',psptyp,'Z',Z,'mass',mass,'coords',coords,'mag',mag,'n_atm_typ',n_atm_typ,'lloc',lloc);
end
fclose(fid1);
assert(size(Atoms,1)==n_atm,'Error: Number of atoms read does not match that specified!');

fprintf(' Finished reading input file!\n');
% disp('###################################################################');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Initialization                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Starting initialization of quantities...\n');

option = C_meth{1,1};

% Lattice unit vectors and Jacobian
lat_vec = zeros(3,3);
lat_vec(1,:) = cell2mat(C_vec1);
lat_vec(2,:) = cell2mat(C_vec2);
lat_vec(3,:) = cell2mat(C_vec3);

Jacb = (det(lat_vec'));
assert(Jacb > 0.0,'Volume is negative!')

% metric_T, Gradient and laplacian transformation matrices
metric_T = lat_vec * lat_vec' ;
metric_T(1,2) = 2*metric_T(1,2); metric_T(2,3) = 2*metric_T(2,3); metric_T(1,3) = 2*metric_T(1,3);
grad_T = inv(lat_vec') ;
lapc_T = grad_T * grad_T' ;
lapc_T(1,2) = 2*lapc_T(1,2); lapc_T(2,3) = 2*lapc_T(2,3); lapc_T(1,3) = 2*lapc_T(1,3);
%disp(lapc_T);aa
% Simulation domain
L1 = sides{1,1};
L2 = sides{1,2};
L3 = sides{1,3};

% Boundary condition type: 1--isolated cluster; 2--periodic system
BC = C_BC{1,1};
if (option == 2) && ((BC == 1) || (BC == 4))
    error('Non-orthogonal option allowed only for boundary conditions 2 & 3');
end

% Finite-difference discretization
if(BC == 1)
    Nx = C_int{1,1} + 1;
    Ny = C_int{1,2} + 1;
    Nz = C_int{1,3} + 1;
elseif(BC == 2)
    Nx = C_int{1,1};
    Ny = C_int{1,2};
    Nz = C_int{1,3};
elseif(BC == 3)
    Nx = C_int{1,1};
    Ny = C_int{1,2};
    Nz = C_int{1,3} + 1;
elseif(BC == 4)
    Nx = C_int{1,1} + 1;
    Ny = C_int{1,2} + 1;
    Nz = C_int{1,3};
    
else
    error('Boundary condition should be one among {1,2,3,4}');
end

dx = L1 / C_int{1,1};
dy = L2 / C_int{1,2};
dz = L3 / C_int{1,3};
N = Nx * Ny * Nz; % total number of FD nodes

% Brillouin-Zone Sampling

nkpt = cell2mat(C_kpt);
if(BC == 1)
    if (nkpt(1)>1 || nkpt(2)>1 || nkpt(3) >1)
        error('For boundary condition 1, nkpt cannot be greater than 1');
    end
elseif(BC == 3)
    if (nkpt(3) >1)
        error('For boundary condition 3, nkpt_z cannot be greater than 1');
    end
elseif(BC == 4)
    if (nkpt(1)>1 || nkpt(2)>1)
        error('For boundary condition 4, nkpt_x and nkpt_y cannot be greater than 1');
    end
end
kptgrid_x = (2*(1:1:nkpt(1)) - nkpt(1)*ones(1,nkpt(1)) -ones(1,nkpt(1)))/2/nkpt(1);
kptgrid_y = (2*(1:1:nkpt(2)) - nkpt(2)*ones(1,nkpt(2)) -ones(1,nkpt(2)))/2/nkpt(2);
kptgrid_z = (2*(1:1:nkpt(3)) - nkpt(3)*ones(1,nkpt(3)) -ones(1,nkpt(3)))/2/nkpt(3);
[kptgrid_X, kptgrid_Y, kptgrid_Z] = meshgrid(kptgrid_x,kptgrid_y,kptgrid_z);
kptgrid_X = permute(kptgrid_X,[2 1 3]);
kptgrid_Y = permute(kptgrid_Y,[2 1 3]);
kptgrid_Z = permute(kptgrid_Z,[2 1 3]);
kptgrid = cat(2,reshape(kptgrid_X,[],1),reshape(kptgrid_Y,[],1),reshape(kptgrid_Z,[],1));
tnkpt = prod(nkpt);
wkpt = ones(tnkpt,1)/tnkpt;% weights for k-points

% Time-Reversal Symmetry to reduce k-points
if C_trevsym{1,1} == 1
    Ikpt = zeros(tnkpt,1);
    Ikpt_rev = zeros(tnkpt,1);
    for ii = 1:tnkpt
        for jj = ii+1:tnkpt
            if (abs(kptgrid(ii,1) + kptgrid(jj,1)) < 1e-8) && (abs(kptgrid(ii,2) + kptgrid(jj,2)) < 1e-8) && (abs(kptgrid(ii,3) + kptgrid(jj,3)) < 1e-8)
                Ikpt(ii) = 1;
                Ikpt_rev(jj) = 1;
            end
        end
    end
    Ikpt = Ikpt>0.5;
    Ikpt_rev = Ikpt_rev>0.5;
    wkpt(Ikpt_rev) = 2*wkpt(Ikpt_rev);
    kptgrid = kptgrid(~Ikpt,:);
    wkpt = wkpt(~Ikpt);
    tnkpt = size(wkpt,1);
end

%  tnkpt = 21;
%  kptgrid = [0.00000000E+00  0.00000000E+00  8.33333333E-02
%             1.66666667E-01  0.00000000E+00  8.33333333E-02
%             3.33333333E-01  0.00000000E+00  8.33333333E-02
%             5.00000000E-01  0.00000000E+00  8.33333333E-02
%             1.66666667E-01  1.66666667E-01  8.33333333E-02
%             3.33333333E-01  1.66666667E-01  8.33333333E-02
%             3.33333333E-01  3.33333333E-01  8.33333333E-02
%             0.00000000E+00  0.00000000E+00  2.50000000E-01
%             1.66666667E-01  0.00000000E+00  2.50000000E-01
%             3.33333333E-01  0.00000000E+00  2.50000000E-01
%             5.00000000E-01  0.00000000E+00  2.50000000E-01
%             1.66666667E-01  1.66666667E-01  2.50000000E-01
%             3.33333333E-01  1.66666667E-01  2.50000000E-01
%             3.33333333E-01  3.33333333E-01  2.50000000E-01
%             0.00000000E+00  0.00000000E+00  4.16666667E-01
%             1.66666667E-01  0.00000000E+00  4.16666667E-01
%             3.33333333E-01  0.00000000E+00  4.16666667E-01
%             5.00000000E-01  0.00000000E+00  4.16666667E-01
%             1.66666667E-01  1.66666667E-01  4.16666667E-01
%             3.33333333E-01  1.66666667E-01  4.16666667E-01
%             3.33333333E-01  3.33333333E-01  4.16666667E-01];
%   wkpt = [0.00926    0.05556    0.05556    0.02778    0.05556    0.11111 ...
%           0.01852    0.00926    0.05556    0.05556    0.02778    0.05556 ...
%           0.11111    0.01852    0.00926    0.05556    0.05556    0.02778 ...
%           0.05556    0.11111    0.01852]';
%   tnkpt = 2;
%   kptgrid = [-2.50000000E-01  5.00000000E-01  0.00000000E+00
%              -2.50000000E-01  0.00000000E+00  0.00000000E+00];
%   wkpt = [0.75 0.25];

% tnkpt = 1;
% kptgrid = [0.25 0.25 0.25];
% wkpt = 1;

% Band structure flags
% isBS = C_bs{1,1};
% lattice = C_lattice{1,1};
isBS = 0; % disable band structure
lattice = 'undefined';

% Degree of Chebyshev polynomial
npl = C_npl{1,1};

% Preconditioning of AAJ/PP
conditioner = C_conditioner{1,1};

% Electronic temperature
Temp = C_thet{1,1};

% Maximum number of relaxation steps
max_relax_it = C_maxrlx{1,1};

% Tolerance for SCF convergence
SCF_tol = C_scf{1,1};

% Number of states
Nev=C_nev{1,1};
if Nev < Nelectron/2
    error('No. of states (Nev) must not be less than half of No. of total valence electrons (Nelectron)!');
end

% Tolerance for Poisson solver
poisson_tol = C_poist{1,1};

% Tolerance for calculating rb
pseudocharge_tol = C_pseudchgt{1,1};

% rc reference
rc_ref = C_rcref{1,1};

% Factor for conversion from Ha to eV
Cst = 27.211384523 ;

% Inverse of smearing
bet = Cst/(8.617343e-5*Temp) ;

% Half of finite difference order
FDn = C_fdn{1,1};

% Finite difference weights of the second derivative
w2 = zeros(1,FDn+1) ;
for k=1:FDn
    w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
        (k*k*factorial(FDn-k)*factorial(FDn+k));
    w2(1) = w2(1)-2*(1/(k*k));
end

% Finite difference weights of the first derivative
w1 = zeros(1,FDn) ;
for k=1:FDn
    w1(k+1) = ((-1)^(k+1))*(factorial(FDn)^2)/...
        (k*factorial(FDn-k)*factorial(FDn+k));
end

% Weights for spatial integration over domain
dV = dx * dy * dz*Jacb ;
W = ones(N,1) * dV;
% For zero-dirichlet, boundary nodes have smaller weights
count = 1;
for kk = 1:Nz
    for jj = 1:Ny
        for ii = 1:Nx
            if (ii == 1)
                W(count) = (BC == 1) * 0.5 * W(count) + (BC == 2) * W(count) + (BC == 3) * W(count) + (BC == 4) * 0.5 * W(count);
            end
            if (ii == Nx)
                W(count) = (BC == 1) * 0.5 * W(count) + (BC == 2) * W(count) + (BC == 3) * W(count) + (BC == 4) * 0.5 * W(count);
            end
            if (jj == 1)
                W(count) = (BC == 1) * 0.5 * W(count) + (BC == 2) * W(count) + (BC == 3) * W(count) + (BC == 4) * 0.5 * W(count);
            end
            if (jj == Ny)
                W(count) = (BC == 1) * 0.5 * W(count) + (BC == 2) * W(count) + (BC == 3) * W(count) + (BC == 4) * 0.5 * W(count);
            end
            if (kk == 1)
                W(count) = (BC == 1) * 0.5 * W(count) + (BC == 2) * W(count) + (BC == 3) * 0.5 * W(count) + (BC == 4) * W(count);
            end
            if (kk == Nz)
                W(count) = (BC == 1) * 0.5 * W(count) + (BC == 2) * W(count) + (BC == 3) * 0.5 * W(count) + (BC == 4) * W(count);
            end
            count = count + 1;
        end
    end
end

% Calculate Spherical Harmonics with origin shifted to the center of the domain
xx_aug = [0-FDn:Nx+FDn-1]*dx - L1/2;
yy_aug = [0-FDn:Ny+FDn-1]*dy - L2/2;
zz_aug = [0-FDn:Nz+FDn-1]*dz - L3/2;
[XX_AUG_3D,YY_AUG_3D,ZZ_AUG_3D] = meshgrid(xx_aug,yy_aug,zz_aug);
XX_AUG_3D = permute(XX_AUG_3D,[2 1 3]);
YY_AUG_3D = permute(YY_AUG_3D,[2 1 3]);
ZZ_AUG_3D = permute(ZZ_AUG_3D,[2 1 3]);
if option == 1
    RR_AUG_3D = sqrt(XX_AUG_3D.^2 + YY_AUG_3D.^2 + ZZ_AUG_3D.^2);
else    
    RR_AUG_3D = sqrt(metric_T(1,1)*XX_AUG_3D.^2 + metric_T(1,2)*(XX_AUG_3D.*YY_AUG_3D) + metric_T(1,3)*(XX_AUG_3D.*ZZ_AUG_3D) + ...
                     metric_T(2,2)*YY_AUG_3D.^2 + metric_T(2,3)*(YY_AUG_3D.*ZZ_AUG_3D) + metric_T(3,3)*ZZ_AUG_3D.^2) ;
end
XX_AUG = reshape(XX_AUG_3D,[],1);
YY_AUG = reshape(YY_AUG_3D,[],1);
ZZ_AUG = reshape(ZZ_AUG_3D,[],1);
RR_AUG = reshape(RR_AUG_3D,[],1);


% isIn = (XX_AUG >= -L1/2) & (XX_AUG <= L1/2) & (YY_AUG >= -L2/2) &...
in_flag = ones(Nx+2*FDn,Ny+2*FDn,Nz+2*FDn);
in_flag(1:FDn,:,:) = 0;
in_flag(:,1:FDn,:) = 0;
in_flag(:,:,1:FDn) = 0;
in_flag(FDn+Nx+1:end,:,:) = 0;
in_flag(:,FDn+Ny+1:end,:) = 0;
in_flag(:,:,FDn+Nz+1:end) = 0;
isIn = (in_flag~=0);

% isIn_3D = reshape(isIn,size(RR_AUG_3D));
% isIn_3D = (XX_AUG_3D >= -L1/2) & (XX_AUG_3D <= L1/2) & (YY_AUG_3D >= -L2/2) &...
% (YY_AUG_3D <= L2/2) & (ZZ_AUG_3D >= -L3/2) & (ZZ_AUG_3D <= L3/2);

RR = RR_AUG(isIn);
% xx = [0:Nx-1]*dx - L1/2;
% yy = [0:Ny-1]*dy - L2/2;
% zz = [0:Nz-1]*dz - L3/2;
% [XX_3D,YY_3D,ZZ_3D] = meshgrid(xx,yy,zz);
% XX_3D = permute(XX_3D,[2 1 3]);
% YY_3D = permute(YY_3D,[2 1 3]);
% ZZ_3D = permute(ZZ_3D,[2 1 3]);
% RR_3D = sqrt(XX_3D.^2 + YY_3D.^2 + ZZ_3D.^2);
% XX = reshape(XX_3D,[],1);
% YY = reshape(YY_3D,[],1);
% ZZ = reshape(ZZ_3D,[],1);
% RR = reshape(RR_3D,[],1);


l_cut = 6;
for l = 0:l_cut
    for m = -l:l
        SH(l+1).Ylm_AUG(:,m+l+1) = SphericalHarmonics(lat_vec,XX_AUG,YY_AUG,ZZ_AUG,l,m,'real');
        Ylm_AUG_TEMP = SH(l+1).Ylm_AUG(:,m+l+1);
        SH(l+1).Ylm(:,m+l+1) = Ylm_AUG_TEMP(isIn);
        % SH(l+1).Ylm(:,m+l+1) = SphericalHarmonics(XX,YY,ZZ,l,m,'real');
    end
end


% Creating structure of all the variables
S = struct('option',option,'lat_vec',lat_vec,'metric_T',metric_T,'grad_T',grad_T,'lapc_T',lapc_T,'Jacb',Jacb,'L1',L1,'L2',L2,'L3',L3,'Nx',Nx,'Ny',Ny,'Nz',Nz,'kptgrid',kptgrid,...
    'tnkpt',tnkpt,'wkpt',wkpt,'BC',BC,'isBS',isBS,'lattice',lattice,...
    'dx',dx,'dy',dy,'dz',dz,'SCF_tol',SCF_tol,'Nev',Nev,...
    'poisson_tol',poisson_tol,'pseudocharge_tol',pseudocharge_tol,...
    'Atoms',Atoms,'n_atm',n_atm,...
    'Temp',Temp,'bet',bet,'npl',npl,'conditioner',conditioner,'FDn',FDn,'rc_ref',rc_ref,...
    'w2',w2,'W',W,'Cst',Cst,'max_relax_it',max_relax_it,...
    'w1',w1,'dV',dV,'dbg_switch',dbg_switch,'nspin', nspin,'xc',xc,...
    'N',N,'Nelectron',Nelectron,'n_typ',n_typ,'Atm',Atm,'fname', filename,...
    'isIn',isIn,'l_cut',l_cut,'RR',RR,'RR_AUG_3D',RR_AUG_3D,'SH',SH);
% Save the initial positions, since later S.Atoms will be changed
S.Atoms_init = S.Atoms;

% Store the factor with which occupation numbers will be multiplied
S.occfac = (3 - S.nspin); % 1 = spin polarized and 2 - spin unpolarized
fprintf(' Finished initialization of quantities!\n');
% disp('###################################################################');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************************************************************
%*                      Make Pseudopotential files                        *
%**************************************************************************
% generate the Pseudopotential_file.mat file for the current atom and Vloc
fprintf(' Reading Pseudopotential files ...\n');

dx2 = S.dx*S.dx;
dy2 = S.dy*S.dy;
dz2 = S.dz*S.dz;
dxdy = S.dx*S.dy;
dydz = S.dy*S.dz;
dzdx = S.dz*S.dx;

% Starting and ending indices of b-region
rb_up = 15;
ii_s_temp = -ceil(rb_up/S.dx);
ii_e_temp = ceil(rb_up/S.dx);
jj_s_temp = -ceil(rb_up/S.dy);
jj_e_temp = ceil(rb_up/S.dy);
%kk_s_temp = -ceil(rb_up/S.dz);
kk_s_temp = 0;
kk_e_temp = ceil(rb_up/S.dz);
xx_temp = (ii_s_temp-S.FDn:ii_e_temp+S.FDn)*S.dx;
yy_temp = (jj_s_temp-S.FDn:jj_e_temp+S.FDn)*S.dy;
zz_temp = (kk_s_temp-S.FDn:kk_e_temp+S.FDn)*S.dz;
[XX_3D_temp,YY_3D_temp,ZZ_3D_temp] = meshgrid(xx_temp,yy_temp,zz_temp);
XX_3D_temp = permute(XX_3D_temp,[2 1 3]);
YY_3D_temp = permute(YY_3D_temp,[2 1 3]);
ZZ_3D_temp = permute(ZZ_3D_temp,[2 1 3]);
W_temp = ones(size(XX_3D_temp)) * S.dV ;
xx_s = ii_s_temp*S.dx; yy_s = jj_s_temp*S.dy; zz_s = kk_s_temp*S.dz;
xx_e = ii_e_temp*S.dx; yy_e = jj_e_temp*S.dy; zz_e = kk_e_temp*S.dz;
W_temp(XX_3D_temp == xx_s) = 0.5 * W_temp(XX_3D_temp == xx_s);
W_temp(YY_3D_temp == yy_s) = 0.5 * W_temp(YY_3D_temp == yy_s);
W_temp(ZZ_3D_temp == zz_s) = 0.5 * W_temp(ZZ_3D_temp == zz_s);
W_temp(XX_3D_temp == xx_e) = 0.5 * W_temp(XX_3D_temp == xx_e);
W_temp(YY_3D_temp == yy_e) = 0.5 * W_temp(YY_3D_temp == yy_e);
W_temp(ZZ_3D_temp == zz_e) = 0.5 * W_temp(ZZ_3D_temp == zz_e);
if S.option == 1
    dd_temp = sqrt(XX_3D_temp.^2 + YY_3D_temp.^2 + ZZ_3D_temp.^2);
else
    dd_temp = sqrt(S.metric_T(1,1)*XX_3D_temp.^2 + S.metric_T(1,2)*(XX_3D_temp.*YY_3D_temp) + S.metric_T(1,3)*(XX_3D_temp.*ZZ_3D_temp) + ...
                   S.metric_T(2,2)*YY_3D_temp.^2  + S.metric_T(2,3)*(YY_3D_temp.*ZZ_3D_temp) + S.metric_T(3,3)*ZZ_3D_temp.^2);
end
for ii = 1:n_typ
    element = Atm(ii).typ;
    lloc = Atm(ii).lloc;
    if Atm(ii).psptyp == 0
        ReadPseudoPot(lloc,element);
    elseif Atm(ii).psptyp == 1
        ReadPseudo_ONCV_Pot(element);
    end
    
    %**********************************************************************
    %*                            Calculate rb                            *
    %**********************************************************************
    if S.Atm(ii).psptyp == 0
        filename = strcat('./Pseudopotentials/Pseudopotential_', S.Atm(ii).typ);
        load(filename)
    elseif S.Atm(ii).psptyp == 1
        filename = strcat('./Pseudopotentials/Pseudopotential_oncv_', S.Atm(ii).typ);
        load(filename)
    end
    %filename = strcat('./Pseudopotentials/Pseudopotential_', S.Atm(ii).typ);
    %load(filename)
    %V_PS_temp = interp1(r_grid_vloc, r_grid_vloc.*Vloc, dd_temp, 'spline');
    V_PS_temp = zeros(size(dd_temp));
    IsLargeThanRmax = dd_temp > r_grid_vloc(end);
    V_PS_temp(IsLargeThanRmax) = -S.Atm(ii).Z;
    V_PS_temp(~IsLargeThanRmax) = interp1(r_grid_vloc, r_grid_vloc.*Vloc, dd_temp(~IsLargeThanRmax), 'spline');
    
    V_PS_temp = V_PS_temp./dd_temp;
    V_PS_temp(dd_temp<r_grid_vloc(2)) = Vloc(1);
    II_temp = 1+S.FDn : size(V_PS_temp,1)-S.FDn;
    JJ_temp = 1+S.FDn : size(V_PS_temp,2)-S.FDn;
    KK_temp = 1+S.FDn : size(V_PS_temp,3)-S.FDn;
    b_temp = zeros(size(V_PS_temp));
    if S.option == 1
        b_temp(II_temp,JJ_temp,KK_temp) = S.w2(1) * (1/dx2 + 1/dy2 + 1/dz2) * V_PS_temp(II_temp,JJ_temp,KK_temp);
        for p = 1:S.FDn
            b_temp(II_temp,JJ_temp,KK_temp) = b_temp(II_temp,JJ_temp,KK_temp) + S.w2(p+1)/dx2 * (V_PS_temp(II_temp+p,JJ_temp,KK_temp) + V_PS_temp(II_temp-p,JJ_temp,KK_temp)) + ...
                                            + S.w2(p+1)/dy2 * (V_PS_temp(II_temp,JJ_temp+p,KK_temp) + V_PS_temp(II_temp,JJ_temp-p,KK_temp)) + ...
                                            + S.w2(p+1)/dz2 * (V_PS_temp(II_temp,JJ_temp,KK_temp+p) + V_PS_temp(II_temp,JJ_temp,KK_temp-p));
        end
    else
        b_temp(II_temp,JJ_temp,KK_temp) = S.w2(1) * (S.lapc_T(1,1)/dx2 + S.lapc_T(2,2)/dy2 + S.lapc_T(3,3)/dz2) * V_PS_temp(II_temp,JJ_temp,KK_temp);
        for p = 1:S.FDn
            b_temp(II_temp,JJ_temp,KK_temp) = b_temp(II_temp,JJ_temp,KK_temp) + S.w2(p+1)*S.lapc_T(1,1)/dx2 * (V_PS_temp(II_temp+p,JJ_temp,KK_temp) + V_PS_temp(II_temp-p,JJ_temp,KK_temp)) + ...
                + S.w2(p+1)*S.lapc_T(2,2)/dy2 * (V_PS_temp(II_temp,JJ_temp+p,KK_temp) + V_PS_temp(II_temp,JJ_temp-p,KK_temp)) + ...
                + S.w2(p+1)*S.lapc_T(3,3)/dz2 * (V_PS_temp(II_temp,JJ_temp,KK_temp+p) + V_PS_temp(II_temp,JJ_temp,KK_temp-p));
            for q = 1:S.FDn
                b_temp(II_temp,JJ_temp,KK_temp) = b_temp(II_temp,JJ_temp,KK_temp) + S.w1(p+1)*S.w1(q+1)*S.lapc_T(1,2)/dxdy * ( V_PS_temp(II_temp+q,JJ_temp+p,KK_temp) ...
                    - V_PS_temp(II_temp-q,JJ_temp+p,KK_temp) - V_PS_temp(II_temp+q,JJ_temp-p,KK_temp) + V_PS_temp(II_temp-q,JJ_temp-p,KK_temp) ) ...
                    + S.w1(p+1)*S.w1(q+1)*S.lapc_T(2,3)/dydz * ( V_PS_temp(II_temp,JJ_temp+q,KK_temp+p) ...
                    - V_PS_temp(II_temp,JJ_temp-q,KK_temp+p) - V_PS_temp(II_temp,JJ_temp+q,KK_temp-p) + V_PS_temp(II_temp,JJ_temp-q,KK_temp-p) ) ...
                    + S.w1(p+1)*S.w1(q+1)*S.lapc_T(1,3)/dzdx * ( V_PS_temp(II_temp+q,JJ_temp,KK_temp+p) ...
                    - V_PS_temp(II_temp-q,JJ_temp,KK_temp+p) - V_PS_temp(II_temp+q,JJ_temp,KK_temp-p) + V_PS_temp(II_temp-q,JJ_temp,KK_temp-p) ) ;
            end
        end
    end
    b_temp = -b_temp / (4*pi);
    err_rb = 100;
    count = 1;
    rb = rc;
    while (err_rb > S.pseudocharge_tol && count <= 100 && rb <= rb_up)
        ii_rb = -1*ii_s_temp+S.FDn-ceil(rb/S.dx)+1:-1*ii_s_temp+S.FDn+ceil(rb/S.dx)+1;
        jj_rb = -1*jj_s_temp+S.FDn-ceil(rb/S.dy)+1:-1*jj_s_temp+S.FDn+ceil(rb/S.dy)+1;
        kk_rb = S.FDn+1:S.FDn+ceil(rb/S.dz)+1;
        rb = rb + max(max(S.dx,S.dy),S.dz);
        err_rb = abs(sum(sum(sum(W_temp(ii_rb,jj_rb,kk_rb).*b_temp(ii_rb,jj_rb,kk_rb)))) + S.Atm(ii).Z/2)/(S.Atm(ii).Z/2);
        count = count + 1;
    end
    %S.lapc_T
    %S.Atm(ii).Z
    assert(rb<=rb_up,'Need to increase upper bound for rb!');
    rb
    %S.Atm(ii).rb = 7;
    S.Atm(ii).rb = rb;
    %S.Atm(ii).rb
    %S.Atm(ii).rb = 9;
    %S.Atm(1).rb = 4.73;
    %S.Atm(2).rb = 4.33;
    %	fprintf(2,' rb for atom typ %s is: %f\n',S.Atm(ii).typ,rb)
    
    
    S.Atm(ii).rc = rc;
end
fprintf(' Finished reading pseudopotential files!\n');

end
