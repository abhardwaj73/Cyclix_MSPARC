%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%% INPUTS TO THE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear
%natom = 8;
natom = length(xnew(:,1));
gorder = 100;
vacuum = 11;
% natom is the number of atom in the unit rthogonal cell anf gorder is the
% cyclic group order
% x_orth = [6.5401    1.7832   22.8352
%     1.3625    5.3495   17.2297
%     1.3630    1.7832   21.0022
%     6.5395    5.3495   19.0627
%     3.6063    1.7832   14.1527
%     4.2962    5.3495   25.9122
%    -0.8839    1.7832   14.1556
%     8.7865    5.3495   25.9093]/1.88973;
x_orth = xnew;
x_orth = x_orth';
%x_orth = [0 1.42*cos(pi/3) 1.42+1.42*cos(pi/3) 2*1.42; 0 1.42*sin(pi/3) 1.42*sin(pi/3) 0; 0 0 0 0];%rand(3,natom);
%lat_cons = [3*1.42, sqrt(3)*1.42];%rand(1,2);
%lat_cons = [5.47919294803771 3.77443917930273];
lat_cons = latcons;
% stores the positions of atoms in orthogonal cell
% first row is x-direction and second row is y-direction and third row is
% the out of plane coordinate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcyl1 = zeros(3,natom);
xcyl2 = zeros(3,natom);
% xcyl stores in (r,\theta, z) format
% uses arc assumption to map the coordinates
for i = 1:natom
    x_coord = x_orth(1,i);
    y_coord = x_orth(2,i);
    z_coord = x_orth(3,i);
    % in x-direction
    %xcyl1(1,i) = gorder*lat_cons(1,1)/(2*pi)+z_coord;
    xcyl1(1,i) = lat_cons(1,1)/(2*sin(pi/gorder))+z_coord;
    xcyl1(2,i) = (2*pi/gorder) * (x_coord/lat_cons(1,1));
    xcyl1(3,i) = y_coord;
    %  in y-direction
    %xcyl2(1,i) = gorder*lat_cons(1,2)/(2*pi)+z_coord;
    xcyl2(1,i) = lat_cons(1,2)/(2*sin(pi/gorder))+z_coord;
    xcyl2(2,i) = (2*pi/gorder) * (y_coord/lat_cons(1,2));
    xcyl2(3,i) = x_coord;
end
z_period1 = lat_cons(1,2);
z_period2 = lat_cons(1,1);
% Now generating positions for ABINIT

xabinit1 = zeros(3,gorder*natom);
xabinit2 = zeros(3,gorder*natom);

for i = 1:gorder
    for j = 1:natom
        % in x-direction
         xabinit1(1,(i-1)*natom+j) = xcyl1(1,j)*cos(xcyl1(2,j)+(i-1)*2*pi/gorder);
         xabinit1(2,(i-1)*natom+j) = xcyl1(1,j)*sin(xcyl1(2,j)+(i-1)*2*pi/gorder);
         xabinit1(3,(i-1)*natom+j) = xcyl1(3,j);
       % in y-direction
         xabinit2(1,(i-1)*natom+j) = xcyl2(1,j)*cos(xcyl2(2,j)+(i-1)*2*pi/gorder);
         xabinit2(2,(i-1)*natom+j) = xcyl2(1,j)*sin(xcyl2(2,j)+(i-1)*2*pi/gorder);
         xabinit2(3,(i-1)*natom+j) = xcyl2(3,j);
    end
end

R1 = 2*vacuum + max(xcyl1(1,:)) - min(xcyl1(1,:));
R2 = 2*vacuum + max(xcyl2(1,:)) - min(xcyl2(1,:));
theta1 = 2*pi/gorder;
theta2 = 2*pi/gorder;
z1 = lat_cons(1,2);
z2 = lat_cons(1,1);

data1.X = (xabinit2(1,:)); data1.Y = (xabinit2(2,:)); data1.Z = (xabinit2(3,:)); % For direction1
% data1.X = (xabinit2(1,:)); data1.Y = (xabinit2(2,:)); data1.Z =
% (xabinit2(3,:));  For direction2
mat2pdb(data1)
%molviewer('mat2PDB.pdb')
         
    


xcyl2
cell
mesh
