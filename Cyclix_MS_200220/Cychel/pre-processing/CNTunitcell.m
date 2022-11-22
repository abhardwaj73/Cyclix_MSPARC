% Script to generate atom positions for:
% (1) Unit cell of nanotube of (m,n) chirality
% (2) Unit cell for the bent tube with (m,n) chirality

% Written by Shashikant Kumar, PhD Georgia Tech, 2019
 % Not verified for chiral nanotube (m not equal n)
clc
clear

% inputs
% (1) m, n         --> chiral indices of the nanotube
% (2) N_repeat     --> Number of unit cell in the bent tube, this value
%                      determines the radius of curvature of the bent tube
% (3) sizetori     --> number of unit cell to be plotted of the tori, for
%                      the complete circle, sizetori = N_repeat-1, and for only one cell,
% sizetori = 0;
% (4) bond:        --> Bond length of C-C in 2D graphene sheet
% (5) N;           --> searches for atoms in unit cell from (-N*a1-N*a2) to
%                      (N*a1+N*a2) If some atoms go missing increase N.

%% Input to the code
m =100; n= 9;
N_repeat =30;
sizetori =0;
bond = 2.66; % in Bohr
%bond = .1418;
%bond = 2.660269226782926*0.529177;
%bond = 2.660269226782926+0.05*2.660269226782926;[xcyl,xabinit,xtorus] = cntgen(2.6*0.52,10,0,50,1)
N = 100;
%% Basic quantities calculated 

% lattice vectors
a1 = bond * [sqrt(3);0];
a2 = -bond/2 * [-sqrt(3);3];
atom2 = bond *[0;1];
%chiral vector
chvec = m*a1+n*a2;
% angle of chiral vector and x-axis
angle = atan(chvec(2)/chvec(1));


%diameter of the tube angstrom
dt = sqrt(3)*bond/pi * sqrt(n^2+m^2+n*m)
% chiral angle between 0 and 30 degrees
theta =  atan(sqrt(3) * m/(2*n+m)) * 180/pi
% Translation vectors T = t1a1+t2a2
 t2 = -(2*m+n)/gcd(2*m+n,2*n+m);
 t1 = (2*n+m)/gcd(2*m+n,2*n+m);
% Magniude of T
T = sqrt(3)*bond*sqrt(3*n^2+3*m*n+3*m^2)/gcd(2*m+n,2*n+m)
% Size of unit cell, N is the number of hexagons
Nhexagon = 2*(n^2+m^2+m*n)/gcd(2*m+n,2*n+m);



% twist



alpha_tw = 2*gcd(2*m+n,2*n+m)/sqrt(3)/dt
ncells = T/3/bond
%% Unit cell in 2D sheet for (m,n) tube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Unit cell in 2D Graphene sheet %%%%%%%%%%%%%%%%%%%%%%
sizecell = 2*Nhexagon;
  count = 0;
  count1 = 0;
coord2D = zeros(sizecell,2);
for i = -N:N
    for j = -N:N
       ch = m*a1+n*a2;
       t = t1*a1+t2*a2;
       x = i*a1+j*a2;
     
        if x'*ch>1e-10 && x'*ch <ch'*ch-0.000001 && x'*t>1e-10 && x'*t<t'*t-0.0000001
            count = count+1;
            count1 = count1+1;
            coord_temp = i*a1+j*a2;
            coord2D(count,1) = coord_temp(1);
            coord2D(count,2) = coord_temp(2);
            
        end        
    end
end

count2 = 0;
for i = -N:N
    for j = -N:N
        ch = m*a1+n*a2;
        t = t1*a1+t2*a2;
        x =  atom2+i*a1+j*a2;
        if x'*ch>1e-10 && x'*ch <ch'*ch-0.00001 && x'*t>1e-10 && x'*t<t'*t-0.00001
            count = count+1;
            count2 = count2+1;
            coord_temp = i*a1+j*a2+atom2;
            coord2D(count,1) = coord_temp(1);
            coord2D(count,2) = coord_temp(2);
            
    
        end        
    end
end

count3 = 0;
for i = -N:N
    for j = -N:N
        ch = m*a1+n*a2;
        t = t1*a1+t2*a2;
        x = i*a1+j*a2;
        if x'*ch>-1e-10 && x'*ch <ch'*ch+0.00001 && x'*t<1e-10 && x'*t>-1e-10
            if x'*ch < ch'*ch+0.00001 && x'*ch > ch'*ch-0.00001
            else
                count = count+1;
                count3 = count3+1;            
            coord_temp = i*a1+j*a2;
            coord2D(count,1) = coord_temp(1);
            coord2D(count,2) = coord_temp(2);
            
            end
        end        
    end
end
count4 = 0;
for i = -N:N
    for j = -N:N
        ch = m*a1+n*a2;
        t = t1*a1+t2*a2;
        x =  atom2+i*a1+j*a2;
        if x'*ch>-1e-10 && x'*ch <ch'*ch+0.00001 && x'*t<1e-10 && x'*t>-1e-10
            if x'*ch < ch'*ch+0.00001 && x'*ch > ch'*ch-0.00001
            else
                 count = count+1;
                count4 = count4+1;          
            coord_temp = i*a1+j*a2+atom2;
            coord2D(count,1) = coord_temp(1);
            coord2D(count,2) = coord_temp(2);
           
            end
        end        
    end
end

count5 = 0;
for i = -N:N
    for j = -N:N
        ch = m*a1+n*a2;
        t = t1*a1+t2*a2;
        x = i*a1+j*a2;
        if x'*ch>-1e-10 && x'*ch <1e-10 && x'*t<t'*t+0.0001 && x'*t>-1e-10
            if x'*t < t'*t+0.00001 && x'*t> t'*t-0.00001
            elseif x'*t>-0.0001 && x'*t <0.0001
            else
            count = count+1;
            count5 = count5+1;         
            coord_temp = i*a1+j*a2;
            coord2D(count,1) = coord_temp(1);
            coord2D(count,2) = coord_temp(2);
            
            end
        end        
    end
end

count6 = 0;
for i = -N:N
    for j = -N:N
        ch = m*a1+n*a2;
        t = t1*a1+t2*a2;
        x =  atom2+i*a1+j*a2;
        if x'*ch>-1e-10 && x'*ch <1e-10 && x'*t<t'*t+0.0001 && x'*t>-1e-10
            if x'*t < t'*t+0.00001 && x'*t> t'*t-0.00001[xcyl,xabinit,xtorus] = cntgen(2.6*0.52,10,0,50,1)
            else
                count = count+1;
                count6 = count6+1;            
            coord_temp = i*a1+j*a2+atom2;
            coord2D(count,1) = coord_temp(1);
            coord2D(count,2) = coord_temp(2);
            
            end
        end        
    end
end
%%%%%%%%%%%%%%%%%%%%% Unit cell in 2D graphene sheet %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unit cell of the nanotube in 3D ( wrapping the 2D sheet calculated in last section)

% Transformation of coordinates in 2D cell along chiral and translation
% vectors
transf_mat = [cos(angle) sin(angle); -sin(angle) cos(angle)];
coord2Dcart = coord2D;

coord2D = (transf_mat*coord2D')';
coord2D = fliplr(coord2D);
coordbent = zeros(sizecell,3);
coordCNT = zeros(sizecell,3);

% Wrapping of the transformed coordinates
for i = 1:sizecell
    xtemp = coord2D(i,1);
    ytemp = coord2D(i,2);
    coordCNT(i,1) = xtemp;
    coordCNT(i,2) =  dt/2 * sin(ytemp*2/dt);
    coordCNT(i,3) = dt/2 * cos(2*ytemp/dt);
end

% Superimposing the unit cell to generate the longer tube
if N_repeat>=1
tmod = sqrt(t'*t);
ptemp = coordCNT;
pstore = ptemp;
for i = 1:sizetori
	ptemp = ptemp+tmod*[ones(length(coord2D),1) zeros(length(coord2D),1) zeros(length(coord2D),1)];
	pstore = [pstore;ptemp];
end
else
    pstore = coordCNT;
end
%% Bending of the nanotube

% Radius of the bent tube (depends of the N_repeat)
r_ring = N_repeat*norm(t)/2/pi;
% curvature
k = 1/r_ring;

% bending of the tube unitcell
for i = 1:sizecell
      xtemp = coordCNT(i,1);
      ytemp = coordCNT(i,2);
      ztemp = coordCNT(i,3);
      Rtemp = (1/k)-ytemp;
      scale = 2*pi*Rtemp/(2*pi*(1/k));
      coordbent(i,3) = ztemp;
      coordbent(i,1) = Rtemp*sin(xtemp*scale/Rtemp);
      coordbent(i,2) = Rtemp*cos(xtemp*scale/Rtemp);

end

% bending of the longer tube
bentplot = zeros(length(pstore),3);

for i = 1:length(pstore)
      xtemp = pstore(i,1);
      ytemp = pstore(i,2);
      ztemp = pstore(i,3);
      Rtemp = (1/k)-ytemp;
      scale = 2*pi*Rtemp/(2*pi*(1/k));
      % scale ensures that the tube completes because otherwise outer part
      % of the tube will fall short and inner part will overlap and the
      % tube will not complete.
      bentplot(i,3) = ztemp;
      bentplot(i,1) = Rtemp*sin(xtemp*scale/Rtemp);
      bentplot(i,2) = Rtemp*cos(xtemp*scale/Rtemp);

end

% bent tube in cylindrical coordinates
cylcoord = zeros(length(pstore),3);
for i = 1: length(pstore)
    cylcoord(i,1) = sqrt((bentplot(i,1))^2 + (bentplot(i,2))^2);
    cylcoord(i,2) = atan((bentplot(i,1))/(bentplot(i,2))) * 180/pi;
    cylcoord(i,3) = bentplot(i,3);
end
mn = min(cylcoord(:,3));
cylcoord(:,3) = cylcoord(:,3) - mn + 0.001;

cin = [(1:length(bentplot))' ones(length(bentplot),1) bentplot];
%% Plots

% plot of unit cell of nanotube in 2D sheet
% figure(1)
% t11 = coord2Dcart(:,1);
% t22 = coord2Dcart(:,2);
% plot(t11,t22,'*')
% axis equal


% genetaing the PDB file to print bent nanotube or straight nanotube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%  Printing the PDB file for CNT %%%%%%%%%%%%%%%%%%%%%%%%
% for bent tube
% data1.X = (bentplot(:,1))'; data1.Y = (bentplot(:,2))'; data1.Z = (bentplot(:,3))'; 
% % for straight tube
%   data1.X = (pstore(:,1))'; data1.Y = (pstore(:,2))'; data1.Z = (pstore(:,3))'; 
%  % data1.X = (coordCNT(:,1))'; data1.Y = (coordCNT(:,2))'; data1.Z = (coordCNT(:,3))'; 
% mat2pdb(data1)
% % % % use molviever to generate the plot
%   molviewer('mat2PDB.pdb')
% clear a1mod a2mod angle angledeg Ch coord_temp count1 count2 count3 count4 count5 count6 i j ptemp ...
%       Rtemp t11 t22 temp x xtemp ytemp ztemp