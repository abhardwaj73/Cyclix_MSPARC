function [p1, p2] = dipolemoment(S)
% Dipole moment
X0 = [0 0 0];
[XX,YY,ZZ] = meshgrid([0:S.Nx-1]*S.dx, [0:S.Ny-1]*S.dy, [0:S.Nz-1]*S.dz);
X = reshape(permute(XX,[2 1 3]), [],1);
Y = reshape(permute(YY,[2 1 3]), [],1);
Z = reshape(permute(ZZ,[2 1 3]), [],1);

R = [X Y Z] -X0;
%p1_x = S.W.*(S.rho+S.b).*R(:,1);
%p1_y = S.W.*(S.rho+S.b).*R(:,2);
%p1_z = S.W.*(S.rho+S.b).*R(:,3);
p1_v = S.W.*(S.rho+S.b).*R;
p1 = norm(p1_v)/ 0.393430307; % in Debye

p2 = sum(S.W.*(S.rho).*R,1);
count_typ = 1;
count_typ_atms = 1;
for JJ_a = 1:S.n_atm % loop over all the atoms
    p2 = p2 - S.Atm(count_typ).Z * (S.Atoms(JJ_a,:)-X0);
    if count_typ_atms == S.Atm(count_typ).n_atm_typ
        count_typ_atms = 1;
        count_typ = count_typ + 1;
    else
        count_typ_atms = count_typ_atms + 1;
    end
end
p2 = norm(p2)/ 0.393430307; % in Debye
end

