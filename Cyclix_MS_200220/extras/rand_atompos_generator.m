% Generate random atom position in Ge Supercell.
clc
clear

% ucell = [10.74 10.74 10.74];
% natom = 8; % No. of atoms in fundamental domain
% atomi = [                0                   0                   0
%         2.685000000000000   2.685000000000000   2.685000000000000
%                         0   5.370000000000000   5.370000000000000
%         2.685000000000000   8.055000000000000   8.055000000000000
%         5.370000000000000                   0   5.370000000000000
%         8.055000000000000   2.685000000000000   8.055000000000000
%         5.370000000000000   5.370000000000000                   0
%         8.055000000000000   8.055000000000000   2.685000000000000];

ucell = [5.47 5.47 8.85];
natom = 2;
atomi = [0 0 0
         3.646666666666667   1.823333333333333   4.425000000000000];
perturb = 10; % maximum amount of perturbation (in percentage of acell)
pert_vec = perturb*ucell/100;
ncell = [4 4 4]; % No. of replica of fundamental domain in each direction
cell = ucell.*ncell;
tnatom = natom*prod(ncell);
atomf = zeros(tnatom,3);

rng('default')
rx = rand(tnatom,1)*2 -1;
ry = rand(tnatom,1)*2 -1;
rz = rand(tnatom,1)*2 -1;

for k=1:ncell(3)
    for j=1:ncell(2)
        for i=1:ncell(1)
            sindex = natom*ncell(2)*ncell(1)*(k-1) + natom*ncell(1)*(j-1) +natom*(i-1) +1;
            eindex = sindex+natom -1 ;
            atomf(sindex:eindex,1) = atomi(:,1) + (i-1)*ucell(1) + rx(sindex:eindex)*pert_vec(1) ;
            atomf(sindex:eindex,2) = atomi(:,2) + (j-1)*ucell(2) + ry(sindex:eindex)*pert_vec(2) ;
            atomf(sindex:eindex,3) = atomi(:,3) + (k-1)*ucell(3) + rz(sindex:eindex)*pert_vec(3) ;
        end
    end
end

atomf(:,1) = mod(atomf(:,1) + cell(1),cell(1));
atomf(:,2) = mod(atomf(:,2) + cell(2),cell(2));
atomf(:,3) = mod(atomf(:,3) + cell(3),cell(3));

atomr = atomf./repmat(cell,tnatom,1);
    