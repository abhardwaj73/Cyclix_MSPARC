function [xcyc,xcyc_red,xcart] = orthNT_generator(cell_x,cell_y,X,n,m,twist_ext,ncell_c,ncell_z,vacc,bond_typ,isPlot)
%%

%==========================================================================
%                         Purpose
%==========================================================================
% Generate atom positions for achiral and chiral nanotubes of 2D systems with orthogonal lattice vectors:
% (1) in helical coord. for the unit cell with cyclic and helical symmetry 
% (2) in cartesian coord. for the unit cell with only periodic symmetry
% along axis
%%

%==========================================================================
%                       Details
%==========================================================================
% (1) Nanotubes are generally identified by a pair of indices (n,m)
% (2) Achiral nanotubes are of two types, namely, a) Armchair b) Zig-zag
% (3) Armchair nanotubes are represented by (0,m)
% (4) Zig-zag nanptubes are represented by (n,0) notation
% (5) Twist (rad/length) can be given for both achiral and chiral tubes
% (6) Screw translation is represented by {psi|tau}
%%

%==========================================================================
%                      Input to the code
%==========================================================================
% (1a) cell_x  :       cell dimensions in x-direction (flat sheet u cell)
%      cell_y  :       cell dimensions in y-direction (flat sheet u cell)
% (1b) X       :       atom positions in the unit cell in cartesian coordinates (nx3)
% (2a) (n,m)   :       indices representing the type of nanotube
% (2b)twist_ext:       gives the applied twist per unit length
% (3) ncell_c  :       number of cells in the cyclic direction (only for
% Cyclix code)
% (4) ncell_z  :       number of cells in the axial direction (only for
% Cyclix code)
% (5) vacc     :       vacuum needed in the radial direction
% (6) bond_typ :       type of bond
%                      1 = linear, 2 = arc
% (7) isPlot   :       plot the nanotube
%%

%==========================================================================
%                      Output from the code
%==========================================================================
% xcyc         :   helical coord. of atoms in ncell_c x ncell_z cells
% xcyc_red     :   reduced coordinates corresponding to xcyc coordinates
% xcart        :   cartesian coord. of atoms in the 1D periodic cell
%%
clc
format long

% Unit conversion
bohr2ang = 0.529177210903;
ang2bohr = 1/bohr2ang;
%%

%==========================================================================
%                        Nanotube characterization
%==========================================================================
d = gcd(n,m);
if cell_x >= cell_y
	c2 = round(cell_x^2/cell_y^2);
	c = sqrt(c2);
	a1 = [c 0]*cell_y; a2 = [0 1]*cell_y;
	dR = d*gcd(m/d,c2);
	t1 = m/dR; t2 = -n*c2/dR;
	N_units = (n^2*c2+m^2)/dR;
	M_1 = (n*c2)*d/m/dR;
else
	c2 = round(cell_y^2/cell_x^2);
	c = sqrt(c2);
	a1 = [1 0]*cell_x; a2 = [0 c]*cell_x;	
	dR = d*gcd(n/d,c2);
	t1 = m*c2/dR; t2 = -n/dR;
	N_units = (m^2*c2+n^2)/dR;
	M_1 = n*d/m/dR;
end

Ch_vec = n*a1+m*a2; 
Ch = norm(Ch_vec);
T_vec = t1*a1+t2*a2;
T = norm(T_vec);

% (1) Diameter of the tube
if (d == 1 && bond_typ == 1)
    error('Bond_typ = 1 not allowed for systems with no cyclic symmetry');
elseif (d < 10 && bond_typ == 1)
    fprintf('Warning: Bond_typ = 1 may cause a lot of deviation from actual shape\n');
end

if bond_typ == 1
    dia = 2*(Ch/2/d)/sin(pi/d);
else
    dia = Ch/pi;
end

fprintf('Radius of the tube is %.15f Bohr\n',dia/2);
fprintf('Angle subtended by the unit cell %.15f\n',2*pi/d);
fprintf('Height of the periodic cell (from helical estimate) is %.15f Bohr\n',T);

% (3) Chiral angle
Ch_theta = acos(dot(a1,Ch_vec)/norm(a1)/Ch);
fprintf('Chiral angle is %.15f radians(%.15f degree)\n',Ch_theta,Ch_theta*180/pi);

% (4) Symmetry operations for 2 atom cells (screw translation)
tau = d*T/N_units;
flag = 0;
if m == 0
    M = 1;
    flag = 1;
else
    for q = 0:m/d
        M = M_1 + N_units*q/m;
        if (round(M) == M)
            flag = 1;
            break;
        end
    end
end
if flag == 1
    psi = 2*pi*M/N_units;
else
    error('Failed to get M\n');
end
fprintf('Screw translation corresponding to unit cell is {%f|%f}\n',psi,tau);

% (5) Twist (per unit length) to be given to the code for running 2-atom cell simulation
twist_int = psi/tau;
fprintf('Twist of the tube corresponding to unit cell %f rad/bohr\n',twist_int);
twist = twist_int+twist_ext;
fprintf('Total twist is %.15f rad/bohr\n',twist);
%%

%==========================================================================
%               Atom positions in rectangular flat sheet
%==========================================================================
cell = [Ch/d tau];
natom_u = size(X,1);
X(:,1) = X(:,1)/cell_x;
X(:,2) = X(:,2)/cell_y;
x_uc = findAtom(X,n,m,a1,a2,t1,t2,Ch,T,cell,N_units);

%%

%==========================================================================
%           Atom positions in cylindrical coordinates (and fractional)
%==========================================================================
rad = dia/2;
if(ncell_c > d)
    error('Number of cells in the cyclic direction (ncell_c) cannot be greater than %d\n',d);
end

xcyc = zeros(natom_u*ncell_c*ncell_z,3);
xcyc(1:natom_u,:) = [rad+x_uc(:,3) x_uc(:,1)*(2*pi/d) x_uc(:,2)*cell(2)];
gen_c = [0 2*pi/d 0];
gen_z = [0 psi tau];
for i = 1:ncell_c-1
    xcyc(natom_u*i+1:natom_u*(i+1),:) = xcyc(1:natom_u,:) + gen_c*i;
end
for i = 1:ncell_z-1
    xcyc(natom_u*ncell_c*i+1:natom_u*ncell_c*(i+1),:) = xcyc(1:natom_u*ncell_c,:) + gen_z*i;
end

% xcyc(:,2) = mod(xcyc(:,2),2*pi);

xmin_at = min(xcyc(:,1));
xmax_at = max(xcyc(:,1));
if(xmin_at - vacc <= 0.1)
    error('vacuum asked too large for the given nanotube!');
end

xcyc_red = xcyc./repmat([xmax_at-xmin_at+2*vacc 2*pi*ncell_c/d cell(2)*ncell_z],size(xcyc,1),1);

fprintf('Fractional coordinates in cylindrical coordinates\n');
disp(xcyc_red);

% Find coordinates corresponding to helical coordinate system
xcyc(:,2) = xcyc(:,2) - twist_int*xcyc(:,3);
xcyc(:,2) = mod(xcyc(:,2),2*pi*ncell_c/d);
xcyc_red(:,2) = xcyc(:,2)/(2*pi*ncell_c/d);
%%

%===============================================================================
%      Atom positions in cartesian coordinates for the 1D periodic cell (only fundamental)        
%===============================================================================
twist = twist_ext+twist_int;

if twist == 0
    ncell_z_cart = 1;
    flag = 1;
else
    flag = 0;
    for i = 1:10000*d
        qq = 2*pi*i/d/twist/tau;
        if abs(round(qq)-qq) < 1e-3
            ncell_z_cart = round(qq);
            flag = 1;
            break;
        end
    end
end
if flag == 1
    fprintf('ncell_z_cart = %d\n',ncell_z_cart);
else
    error('No periodic cell could be found in %f rotation\n',2*pi*10000);
end
xcart = zeros(natom_u*d*ncell_z_cart,3);
for ii = 1:ncell_z_cart
    indx = natom_u*d*(ii-1);
    z0 = (ii-1)*cell(2);
    for i = 1:d
        thetatemp = xcyc(1:natom_u,2) + (i-1)*2*pi/d;
        thetatemp = thetatemp + twist * (z0 + xcyc(1:natom_u,3));        
        xcart(indx+natom_u*(i-1)+1:indx+natom_u*i,1) = xcyc(1:natom_u,1).*cos(thetatemp(:));
        xcart(indx+natom_u*(i-1)+1:indx+natom_u*i,2) = xcyc(1:natom_u,1).*sin(thetatemp(:));
        xcart(indx+natom_u*(i-1)+1:indx+natom_u*i,3) = z0 + xcyc(1:natom_u,3);
    end
end
%%

%=============================================================================================
%                       Input for Cyclix SPARC
%=============================================================================================
fprintf('\n\n********************** SPARC INPUT BEGINS (in helical coordinate system)*************\n\n');
fprintf('CELL: %.15f %.15f %.15f\n',xmax_at-xmin_at+2*vacc,2*pi*ncell_c/d,cell(2)*ncell_z);
fprintf('TWIST_ANGLE: %.15f \n',twist);
fprintf('KPT: %d %d %d\n',1,d,ncell_z_cart);
fprintf('[\b"Warning:Change the kpt in the z-direction appropriately"\n]\b');
fprintf('COORD_FRAC:\n')
disp(xcyc_red);
fprintf('************************** SPARC INPUT ENDS ************************************\n\n');
%%

%========================================================================================
%                        Tube visualization
%========================================================================================
if isPlot == 1
    xcart = xcart*bohr2ang;
    data1.X = (xcart(:,1))'; data1.Y = (xcart(:,2))'; data1.Z = (xcart(:,3))';
    mat2pdb(data1)
    molviewer('mat2PDB.pdb')
    xcart = xcart*ang2bohr;
end

end
%%



function x_uc = findAtom(X,n,m,a1,a2,t1,t2,Ch,T,cell,N_units)
    x_uc = X;
    for atom = 1:size(X,1)
        flag = 0;
        for p = 0:N_units
            for q = -N_units:N_units
                pos_Ch = dot((p+X(atom,1))*a1 + (q+X(atom,2))*a2,n*a1+m*a2)/Ch;
                pos_T = dot((p+X(atom,1))*a1 + (q+X(atom,2))*a2,t1*a1+t2*a2)/T;
                if pos_Ch >= 0 && pos_Ch < cell(1) && pos_T >= 0 && pos_T < cell(2)
                    flag = 1;
                    break;
                end
            end
            if(flag == 1)
                break;
            end
        end
        if flag == 0
            error('Could not find the atom positions in the 2D sheet\n');
        else
            x_uc(atom,1:2) = [pos_Ch/cell(1) pos_T/cell(2)];
        end
    end
end