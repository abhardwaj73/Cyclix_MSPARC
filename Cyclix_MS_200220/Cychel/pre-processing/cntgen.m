function [xcyc,xcyc_red,xcart] = nanotube_generator(a,del,n,m,ncell_c,ncell_z,vacc,bond_typ,cell_typ,isPlot)
%%

%==========================================================================
%                         Purpose
%==========================================================================
% Generate atom positions for achiral nanotubes:
% (1) in cylindrical coord. for the unit cell with cyclic and periodic symmetry 
% (2) in cartesian coord. for the unit cell with only periodic symmetry
% along axis
%%

%==========================================================================
%                       Details
%==========================================================================
% (1) Nanotubes are generally identified by a pair of indices (n,m).
% (2) Achiral nanotubes are of two types, namely, a) Armchair b) Zig-zag
% (3) Armchair nanotubes are represented by (n,n) notation.
% (4) Zig-zag nanptubes are represented by (n,0) notation.
% (5) 4-atom cell typ is available only for achiral tubes.
%%

%==========================================================================
%                      Input to the code
%==========================================================================
% (1a) a       :       X-X bond length
% (1b) del     :       out of plane projection of the atom
% (2) (n,m)    :       indices representing the type of nanotube
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
clear
format long
%%

%==========================================================================
%               Atom positions in rectangular flat sheet
%==========================================================================
if cell_typ == 2
    fprintf('Work in progress!');
else 
    if n == m
        fprintf('Tube is Armchair with cyclic group order %d \n',n);
        cell = [3*a sqrt(3)*a 1];
        x_uc = [0 0 -del/2; 1/3 0 del/2; 1/2 1/2 -del/2; 5/6 1/2 del/2] ;
    elseif n == 0
        fprintf('Tube is Zig-Zag with cyclic group order %d \n',n)
        cell = [sqrt(3)*a 3*a 1];
        x_uc = [0 0 -del/2; 0 1/3 del/2; 1/2 1/2 -del/2; 1/2 5/6 del/2] ;
    else
        error('4-atom unit cell only available for chiral tubes!');
    end
end
%%

%==========================================================================
%                    Nanotube characterization
%==========================================================================
d = gcd(n,m);
dR = gcd(2*n+m, 2*m+n);
A = sqrt(3) * a;
Ch = A * sqrt(n^2 + m^2 + n*m); % Circumference

% (1) Diameter of the tube
if bond_typ == 1
    dia = 2*(Ch/2/d)/sin(pi/d);
else
    dia = Ch/pi;
end

fprintf('Radius of the tube is %15.5f \n',dia/2);
fprintf('Angle subtended by the unit cell %.15f\n',2*pi/d);

% (2) Height of the periodic cell
T = sqrt(3) * Ch /dR;
fprintf('Height of the periodic cell is %15.5f \n',T);

% (3) Chiral angle
Ch_theta = acos((2*n+m)/(2*sqrt(n^2+m^2+n*m)));
fprintf('Chiral angle is %15.5f radians (%15.5f degree)\n',Ch_theta,Ch_theta*180/pi);
%%

%==========================================================================
%           Atom positions in cylindrical coordinates
%==========================================================================
rd = dia/2;
if cell_typ == 2
    fprintf('Work in progress!');
else   
    xcyc = zeros(4*ncell_c*ncell_z,3);
    for i = 1:ncell_z
        for j = 1:ncell_c
            xcyc(1,:) = [rd 0 
        if n == m
            theta = (x_uc/3/a)*(2*pi/m);        % for armchair
        elseif n == 0
            theta = (xatom/sqrt(3)/a)*(2*pi/m);  % for zigzag
        end
        xcyl(i,1) = dt/2;
        xcyl(i,2) = theta;
        xcyl(i,3) = yatom;
    end
end

% Find reduced coordinates for SPARC
xmin_at = min(xcyl(:,1));
xmax_at = max(xcyl(:,1));
if(xmin_at - vac <= 0.1)
    error('vacuum asked too large or radius too small!');
end
range1 = xmax_at - xmin_at + 2*vac;
range2 = 2*pi/m;
range3 = length_tube;
xcyl_red = xcyl./repmat([range1 range2 range3],4,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ncell_z == 1
    twist = 0;
else
    twist = 2*pi/(ncell_z*range3);
end
fprintf('Twist %.15f (rad/Bohr)\n',twist);

% The followng code is written only for (m,0) and (m,m) carbon nanotubes
xabinit = zeros(4*m*ncell_z,3);
%% CNT fundamental cell for abinit, Cartesian coordinates (x, y, z)
for ii = 1:ncell_z
    indx = 4*m*(ii-1);
    z0 = (ii-1) * range3;
    for i = 1:m
        thetatemp = [xcyl(1,2)+(i-1)*2*pi/m; xcyl(2,2)+(i-1)*2*pi/m; xcyl(3,2)+(i-1)*2*pi/m; xcyl(4,2)+(i-1)*2*pi/m;];
        thetatemp = thetatemp + twist * (z0 + xcyl(:,3));
        xabinit(indx+4*(i-1)+1,:) = [xcyl(1,1)*cos(thetatemp(1)); xcyl(1,1)*sin(thetatemp(1)); z0+xcyl(1,3)];
        xabinit(indx+4*(i-1)+2,:) = [xcyl(2,1)*cos(thetatemp(2)); xcyl(2,1)*sin(thetatemp(2)); z0+xcyl(2,3)];
        xabinit(indx+4*(i-1)+3,:) = [xcyl(3,1)*cos(thetatemp(3)); xcyl(3,1)*sin(thetatemp(3)); z0+xcyl(3,3)];
        xabinit(indx+4*(i-1)+4,:) = [xcyl(4,1)*cos(thetatemp(4)); xcyl(4,1)*sin(thetatemp(4)); z0+xcyl(4,3)];
    end
end

%% Fundamental unit cell of torus in cylindrical coordinates (r, theta, z)
% r_torus = length_tube/(2*sin(pi/gorder_torus));
% fprintf('Group order of the torus is %d \n',gorder_torus);
% fprintf('Radius of the torus is %15.5f \n',r_torus);
% xtorus = zeros(4*m,3);
% for i = 1:4*m
%     xtorus(i,1) = r_torus - xabinit(i,1);
%     xtorus(i,2) = xabinit(i,3)/length_tube * 2*pi/gorder_torus * 180/pi;
%     xtorus(i,3) = xabinit(i,2);
% end

% %% coordinates of whole torus in cartesian coordinates (x, y, z)
% xtorus_full = zeros(gorder_torus*4*m,3);

% for i = 1:gorder_torus
%     for j = 1:4*m
%         thetatemp = xtorus(j,2) + (i-1)*2*180/gorder_torus;
%         xtorus_full((i-1)*4*m+j,1) = xtorus(j,1)*cos(thetatemp*pi/180);
%         xtorus_full((i-1)*4*m+j,2) = xtorus(j,1)*sin(thetatemp*pi/180);
%         xtorus_full((i-1)*4*m+j,3) = xtorus(j,3);
%     end
% end

%% Plotting the tube
if ifplot == 1
    data1.X = (xtorus_full(:,1))'; data1.Y = (xtorus_full(:,2))'; data1.Z = (xtorus_full(:,3))';
    %  data1.X = (xbent(:,1))'; data1.Y = (xbent(:,2))'; data1.Z = (xbent(:,3))';
    %data1.X = (xCNT(:,1))'; data1.Y = (xCNT(:,2))'; data1.Z = (xCNT(:,3))';
    mat2pdb(data1)
    molviewer('mat2PDB.pdb')
end
