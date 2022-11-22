function [xcyc,xcyc_red,xcart] = nanotube_generator(a,del,n,m,twist_ext,ncell_c,ncell_z,vacc,bond_typ,cell_typ,isPlot)
%%

%==========================================================================
%                         Purpose
%==========================================================================
% Generate atom positions for achiral and chiral nanotubes:
% (1) in cylindrical coord. for the unit cell with cyclic and periodic symmetry
% (2) in helical coord. for the unit cell with cyclic and helical symmetry 
% (2) in cartesian coord. for the unit cell with only periodic symmetry
% along axis.
%%

%==========================================================================
%                       Details
%==========================================================================
% (1) Nanotubes are generally identified by a pair of indices (n,m).
% (2) Achiral nanotubes are of two types, namely, a) Armchair b) Zig-zag
% (3) Armchair nanotubes are represented by (n,n) notation.
% (4) Zig-zag nanptubes are represented by (n,0) notation.
% (5) 4-atom cell typ is available only for achiral tubes.
% (6) Twist (rad/length) can be given for both achiral and chiral tubes.
% (7) Screw translation is represented by {psi|tau}
%%

%==========================================================================
%                      Input to the code
%==========================================================================
% (1a) a       :       X-X bond length
% (1b) del     :       out of plane projection of the atom
% (2a) (n,m)   :       indices representing the type of nanotube
% (2b)twist_ext:       gives the applied twist per unit length
% (3) ncell_c  :       number of cells in the cyclic direction
% (4) ncell_z  :       number of cells in the axial direction
% (5) vacc     :       vacuum needed in the radial direction
% (6) bond_typ :       type of bond
%                      1 = linear, 2 = arc
% (7) cell_typ :       2 = 2 atom cell, 4 =  4 atom cell
% (8) isPlot   :       plot the nanotube
%%

%==========================================================================
%                      Output from the code
%==========================================================================
% xcyc         :   cylindrical coord. of atoms in ncell_c x ncell_z cells
% xcyc_red     :   reduced coordinates corresponding to xcyc coordinates
% xcart        :   cartesian coord. of atoms
%%
clc
format long

bohr2ang = 0.529177210903;
ang2bohr = 1/bohr2ang;

%==========================================================================
%                        Nanotube characterization
%==========================================================================
A = sqrt(3)*a;
a1 = [sqrt(3)/2 1/2]*A; a2 = [sqrt(3)/2 -1/2]*A;
d = gcd(n,m);
dR = gcd(2*n+m, 2*m+n);
Ch = A * sqrt(n^2 + m^2 + n*m); % Circumference
t1 = (2*m+n)/dR; t2 = -(2*n+m)/dR; 

% (1) Diameter of the tube
if bond_typ == 1
    dia = 2*(Ch/2/d)/sin(pi/d);
else
    dia = Ch/pi;
end

fprintf('Radius of the tube is %.15f Bohr\n',dia/2);
fprintf('Angle subtended by the unit cell %.15f\n',2*pi/d);

% (2) Height of the periodic cell
T = sqrt(3) * Ch /dR;
fprintf('Height of the periodic cell (from helical estimate) is %.15f Bohr\n',T);

% (3) Chiral angle
Ch_theta = acos((2*n+m)/(2*sqrt(n^2+m^2+n*m)));
fprintf('Chiral angle is %.15f radians(%.15f degree)\n',Ch_theta,Ch_theta*180/pi);

% (4) Symmetry operations for 2 atom cells (screw translation)
N_hex = 2*(n^2+m^2+n*m)/dR;
tau = d*T/N_hex;
flag = 0;
if m == 0
    M = 1;
    flag = 1;
else
    M_1 = (2*n+m)*d/m/dR;
    for q = 0:m/d
        M = M_1 + N_hex*q/m;
        if (round(M) == M)
            flag = 1;
            break;
        end
    end
end
if flag == 1
    psi = 2*pi*M/N_hex;
else
    error('Failed to get M\n');
end
fprintf('Screw translation corresponding to 2 atom cell is {%f|%f}\n',psi,tau);

% (5) Twist (per unit length) to be given to the code for running 2-atom cell simulation
twist_int = psi/tau;
fprintf('Twist of the tube corresponding to 2-atom unit cell %f rad/bohr\n',twist_int);
twist = twist_int+twist_ext;
fprintf('Total twist is %.15f rad/bohr\n',twist);
%%

%==========================================================================
%               Atom positions in rectangular flat sheet
%==========================================================================
if cell_typ == 2
    cell = [Ch/d tau];
    %if n == m
    %    pos_Ch = (1/3)*cell(1);
    %    pos_T = 0;
    %else
    %    pos_Ch = ((n+m) + dR*M)*A^2/2/Ch
    %    pos_T = ((m-n)/3 + d)*sqrt(3)*A^2/2/Ch
    %end
    %x_uc = [0 0 -del/2;mod(pos_Ch,cell(1))/cell(1) mod(pos_T,cell(2))/cell(2) del/2];
	
	
	flag = 0;
    for p = 0:N_hex
        for q = -N_hex:N_hex
            pos_Ch = dot((p+1/3)*a1 + (q+1/3)*a2,n*a1+m*a2)/Ch;
            pos_T = dot((p+1/3)*a1 + (q+1/3)*a2,t1*a1+t2*a2)/T;
            if pos_Ch >= 0 && pos_Ch < cell(1) && pos_T >= 0 && pos_T < cell(2)
                flag = 1;
                break;
            end
        end
        if(flag == 1)
            break;
        end
    end
    if(flag == 1)
        x_uc = [0 0 -del/2;pos_Ch/cell(1) pos_T/cell(2) del/2];
    else
        error('Could not find the atom positions in rectangular cell\n');
    end

else 
    if n == m
        cell = [3*a sqrt(3)*a];
        x_uc = [0 0 -del/2; 1/3 0 del/2; 1/2 1/2 -del/2; 5/6 1/2 del/2] ;
    elseif m == 0
        cell = [sqrt(3)*a 3*a 1];
        x_uc = [0 0 -del/2; 0 1/3 del/2; 1/2 1/2 -del/2; 1/2 5/6 del/2] ;
    else
        error('4-atom unit cell only available for achiral tubes of type (n,n) or (n,0)!');
    end
end
%%

%==========================================================================
%           Atom positions in cylindrical coordinates (and fractional)
%==========================================================================
rad = dia/2;
if(ncell_c > d)
    error('Number of cells in the cyclic direction (ncell_c) cannot be greater than %d\n',d);
end
if cell_typ == 2
    xcyc = zeros(2*ncell_c*ncell_z,3);
    xcyc(1:2,:) = [rad+x_uc(:,3) x_uc(:,1)*(2*pi/d) x_uc(:,2)*cell(2)];
    gen_c = [0 2*pi/d 0];
    gen_z = [0 psi tau];
    for i = 1:ncell_c-1
        xcyc(2*i+1:2*(i+1),:) = xcyc(1:2,:) + gen_c*i;
    end
    for i = 1:ncell_z-1
        xcyc(2*ncell_c*i+1:2*ncell_c*(i+1),:) = xcyc(1:2*ncell_c,:) + gen_z*i;
    end
    
    xcyc(:,2) = mod(xcyc(:,2),2*pi);
else   
    xcyc = zeros(4*ncell_c*ncell_z,3);
    xcyc(1:4,:) = [rad+x_uc(:,3) x_uc(:,1)*(2*pi/d) x_uc(:,2)*cell(2)];
    gen_c = [0 2*pi/d 0];
    gen_z = [0 0 cell(2)];
    for i = 1:ncell_c-1
        xcyc(4*i+1:4*(i+1),:) = xcyc(1:4,:) + gen_c*i;
    end

    for i = 1:ncell_z-1
        xcyc(4*ncell_c*i+1:4*ncell_c*(i+1),:) = xcyc(1:4*ncell_c,:) + gen_z*i;
    end
end

xmin_at = min(xcyc(:,1));
xmax_at = max(xcyc(:,1));
if(xmin_at - vacc <= 0.1*bohr2ang)
    error('vacuum asked too large for the given nanotube!');
end

xcyc_red = xcyc./repmat([xmax_at-xmin_at+2*vacc 2*pi*ncell_c/d cell(2)*ncell_z],size(xcyc,1),1);

% Find coordinates corresponding to helical coordinate system
if cell_typ == 2
    xcyc(:,2) = xcyc(:,2) - twist_int*xcyc(:,3);
    xcyc(:,2) = mod(xcyc(:,2),2*pi*ncell_c/d);
    xcyc_red(:,2) = xcyc(:,2)/(2*pi*ncell_c/d);
end
%%

%=============================================================================================
%                       Input for Cyclohelical SPARC
%=============================================================================================
fprintf('\n ***** SPARC INPUT BEGINS ******\n');
fprintf('CELL: %.15f %.15f %.15f\n',xmax_at-xmin_at+2*vacc,2*pi*ncell_c/d,cell(2)*ncell_z);
if cell_typ == 2
    fprintf('TWIST_ANGLE: %.15f \n',twist);
else
    fprintf('TWIST_ANGLE: %.15f \n',twist_ext);
end
fprintf('KPT: %d %d %d\n',1,d,1);
fprintf('Warning:Change the kpt in the z-direction appropriately\n');
fprintf('COORD_FRAC:\n')
disp(xcyc_red);
fprintf('***** SPARC INPUT ENDS ******\n\n');
%%

%===============================================================================
%      Atom positions in cartesian coordinates for the 1D periodic cell (only fundamental)        
%===============================================================================
if cell_typ == 2
    twist = twist_ext+twist_int;
    
    flag = 0;
    for i = 1:10000*d
        qq = 2*pi*i/d/twist/tau;
        if abs(round(qq)-qq) < 1e-3
            ncell_z = round(qq);
            flag = 1;
            break;
        end
    end
    
    if flag == 1
        fprintf('ncell_z = %d\n',ncell_z);
    else
        error('No periodic cell could be found in %f rotation\n',2*pi*10000);
    end
    xcart = zeros(2*d*ncell_z,3);
	for ii = 1:ncell_z
    	indx = 2*d*(ii-1);
    	z0 = (ii-1)*cell(2);
        for i = 1:d
            thetatemp = xcyc(1:2,2) + (i-1)*2*pi/d;
            thetatemp = thetatemp + twist * (z0 + xcyc(1:2,3));
            xcart(indx+2*(i-1)+1,:) = [xcyc(1,1)*cos(thetatemp(1)); xcyc(1,1)*sin(thetatemp(1)); z0+xcyc(1,3)];
            xcart(indx+2*(i-1)+2,:) = [xcyc(2,1)*cos(thetatemp(2)); xcyc(2,1)*sin(thetatemp(2)); z0+xcyc(2,3)];
        end
	end
else
	twist = twist_ext;
    if twist_ext > 0
		%ncell_z = round(2*pi/twist_ext/cell(2));
        for i = 1:d*100
            qq = 2*pi*i/d/twist/cell(2);
            if abs(round(qq)-qq) < 1e-4
                ncell_z = round(qq);
                break;
            end
        end
    fprintf('ncell_z = %d\n',ncell_z);    
    end
    
    xcart = zeros(4*d*ncell_z,3);
	for ii = 1:ncell_z
    	indx = 4*d*(ii-1);
    	z0 = (ii-1)*cell(2);
    	for i = 1:d
        	thetatemp = xcyc(1:4,2) + (i-1)*2*pi/d;
        	thetatemp = thetatemp + twist * (z0 + xcyc(1:4,3));
        	xcart(indx+4*(i-1)+1,:) = [xcyc(1,1)*cos(thetatemp(1)); xcyc(1,1)*sin(thetatemp(1)); z0+xcyc(1,3)];
        	xcart(indx+4*(i-1)+2,:) = [xcyc(2,1)*cos(thetatemp(2)); xcyc(2,1)*sin(thetatemp(2)); z0+xcyc(2,3)];
        	xcart(indx+4*(i-1)+3,:) = [xcyc(3,1)*cos(thetatemp(3)); xcyc(3,1)*sin(thetatemp(3)); z0+xcyc(3,3)];
        	xcart(indx+4*(i-1)+4,:) = [xcyc(4,1)*cos(thetatemp(4)); xcyc(4,1)*sin(thetatemp(4)); z0+xcyc(4,3)];
    	end

	end
end
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
