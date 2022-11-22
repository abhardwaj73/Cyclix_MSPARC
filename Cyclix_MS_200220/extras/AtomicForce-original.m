function force = AtomicForce(S)
% @brief	AtomicForce(S) calculates the atomic force.
%
% @param S				A struct that contains the following (or more) fields:
% @param S.Atm
% @param S.Atoms		A (# of atoms)-by-3 matrix that lists the atom coords of all types
% @param S.dx			The grid size in x-dimesion
% @param S.dy			The grid size in y-dimesion
% @param S.dz			The grid size in z-dimesion
% @param S.FDn			Half of the order of the finite difference scheme
% @param S.Lap			Discretized Laplacian operator (matrix)
% @param S.LapPreconL	Left preconditioner for discretized Laplacian operator (matrix)
% @param S.LapPreconU	Rigtht preconditioner for discretized Laplacian operator (matrix)
% @param S.N 			Total number of grid nodes
% @param S.n_atm		Total number of atoms
% @param S.Nelectron	Number of electrons
% @param S.Nx			The number of nodes in [0,L1]
% @param S.Ny			The number of nodes in [0,L2]
% @param S.Nz			The number of nodes in [0,L3]
% @param S.w1			Finite difference weights of 1st derivative
% @param S.w2			Finite difference weights of 2nd derivative
% @param S.W			Weights for spatial integration over domain
%
% @authors	Qimen Xu <qimenxu@gatech.edu>
%			Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% @2016 Georgia Institute of Technology.
%

fprintf('\n Starting atomic force calculation ... \n');
% tic_forces = tic;

Dpseudo_x = DiscreteGradient(S,0,1)*(S.b + S.b_ref);
Dpseudo_y = DiscreteGradient(S,0,2)*(S.b + S.b_ref);
Dpseudo_z = DiscreteGradient(S,0,3)*(S.b + S.b_ref);

Dphi_x = DiscreteGradient(S,0,1)*(S.phi);
Dphi_y = DiscreteGradient(S,0,2)*(S.phi);
Dphi_z = DiscreteGradient(S,0,3)*(S.phi);
%
% DVc_x = DiscreteGradient(S,0,1)*(S.V_c);
% DVc_y = DiscreteGradient(S,0,2)*(S.V_c);
% DVc_z = DiscreteGradient(S,0,3)*(S.V_c);

% Initialization
%force = zeros(S.n_atm,3);
force_corr = zeros(S.n_atm,3);
force_nloc = zeros(S.n_atm,3);
force_local = zeros(S.n_atm,3);
fac = [S.L1 S.L2 S.L3];
dx2 = S.dx*S.dx;
dy2 = S.dy*S.dy;
dz2 = S.dz*S.dz;
dxdy = S.dx*S.dy;
dydz = S.dy*S.dz;
dzdx = S.dz*S.dx;
coeff = S.w2(1) * (S.lapc_T(1,1)/dx2 + S.lapc_T(2,2)/dy2 + S.lapc_T(3,3)/dz2);

count_typ = 1;
count_typ_atms = 1;
for JJ_a = 1:S.n_atm % loop over all the atoms
    % Load pseudopotential file
    if count_typ_atms == 1
        if S.Atm(count_typ).psptyp == 0
            filename = strcat('./Pseudopotentials/Pseudopotential_', S.Atm(count_typ).typ);
            load(filename)
        elseif S.Atm(count_typ).psptyp == 1
            filename = strcat('./Pseudopotentials/Pseudopotential_oncv_', S.Atm(count_typ).typ);
            load(filename)
        end
        %lloc = S.Atom(JJ_a).lloc;
    end
    % Atom position
    x0 = S.Atoms(JJ_a,1);
    y0 = S.Atoms(JJ_a,2);
    z0 = S.Atoms(JJ_a,3);
    
    if(S.BC == 1)
        % For isolated clusters, each atom has only one image, i.e. itself
        XX_IMG_3D = x0;
        YY_IMG_3D = y0;
        ZZ_IMG_3D = z0;
        n_image_total = 1;
    elseif(S.BC == 2)
        % Note the S.dx, S.dy, S.dz terms are to ensure the image rb-region overlap w/ fund. domain
        n_image_xl = floor((S.Atoms(JJ_a,1) + S.Atm(count_typ).rb)/S.L1);
        n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+S.Atm(count_typ).rb-S.dx)/S.L1);
        n_image_yl = floor((S.Atoms(JJ_a,2) + S.Atm(count_typ).rb)/S.L2);
        n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+S.Atm(count_typ).rb-S.dy)/S.L2);
        n_image_zl = floor((S.Atoms(JJ_a,3) + S.Atm(count_typ).rb)/S.L3);
        n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+S.Atm(count_typ).rb-S.dz)/S.L3);
        % Total No. of images of atom JJ_a (including atom JJ_a)
        n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1) * (n_image_zl+n_image_zr+1);
        % Find the coordinates for all the images
        xx_img = [-n_image_xl : n_image_xr] * S.L1 + x0;
        yy_img = [-n_image_yl : n_image_yr] * S.L2 + y0;
        zz_img = [-n_image_zl : n_image_zr] * S.L3 + z0;
        [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = meshgrid(xx_img,yy_img,zz_img);
        XX_IMG_3D = permute(XX_IMG_3D,[2 1 3]);
        YY_IMG_3D = permute(YY_IMG_3D,[2 1 3]);
        ZZ_IMG_3D = permute(ZZ_IMG_3D,[2 1 3]);
    elseif(S.BC == 3)
        n_image_xl = floor((S.Atoms(JJ_a,1) + S.Atm(count_typ).rb)/S.L1);
        n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+S.Atm(count_typ).rb-S.dx)/S.L1);
        n_image_yl = floor((S.Atoms(JJ_a,2) + S.Atm(count_typ).rb)/S.L2);
        n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+S.Atm(count_typ).rb-S.dy)/S.L2);
        % Total No. of images of atom JJ_a (including atom JJ_a)
        n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1);
        % Find the coordinates for all the images
        xx_img = [-n_image_xl : n_image_xr] * S.L1 + x0;
        yy_img = [-n_image_yl : n_image_yr] * S.L2 + y0;
        zz_img = z0;
        [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = meshgrid(xx_img,yy_img,zz_img);
        XX_IMG_3D = permute(XX_IMG_3D,[2 1 3]);
        YY_IMG_3D = permute(YY_IMG_3D,[2 1 3]);
        ZZ_IMG_3D = permute(ZZ_IMG_3D,[2 1 3]);
    elseif(S.BC == 4)
        n_image_zl = floor((S.Atoms(JJ_a,3) + S.Atm(count_typ).rb)/S.L3);
        n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+S.Atm(count_typ).rb-S.dz)/S.L3);
        
        % Total No. of images of atom JJ_a (including atom JJ_a)
        n_image_total = n_image_zl+n_image_zr+1;
        % Find the coordinates for all the images
        xx_img = x0;
        yy_img = y0;
        zz_img = [-n_image_zl : n_image_zr] * S.L3 + z0;
        [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = meshgrid(xx_img,yy_img,zz_img);
        XX_IMG_3D = permute(XX_IMG_3D,[2 1 3]);
        YY_IMG_3D = permute(YY_IMG_3D,[2 1 3]);
        ZZ_IMG_3D = permute(ZZ_IMG_3D,[2 1 3]);
    end
    
    %     imgRcCount = 0;
    % Loop over all image(s) of atom JJ_a (including atom JJ_a)
    for count_image = 1:n_image_total
        
        % Atom position of the image
        x0_i = XX_IMG_3D(count_image);
        y0_i = YY_IMG_3D(count_image);
        z0_i = ZZ_IMG_3D(count_image);
        
        % Indices of closest grid point to atom
        pos_ii = round(x0_i / S.dx) + 1;
        pos_jj = round(y0_i / S.dy) + 1;
        pos_kk = round(z0_i / S.dz) + 1;
        
        %**********************************************************************
        %*                  Calculate local atomic force term                 *
        %**********************************************************************
        
        % Starting and ending indices of b-region
        ii_s = pos_ii - ceil(S.Atm(count_typ).rb/S.dx+0.5);
        ii_e = pos_ii + ceil(S.Atm(count_typ).rb/S.dx+0.5);
        jj_s = pos_jj - ceil(S.Atm(count_typ).rb/S.dy+0.5);
        jj_e = pos_jj + ceil(S.Atm(count_typ).rb/S.dy+0.5);
        kk_s = pos_kk - ceil(S.Atm(count_typ).rb/S.dz+0.5);
        kk_e = pos_kk + ceil(S.Atm(count_typ).rb/S.dz+0.5);
        if(S.BC == 1)
            % Check if the b-region is inside the domain (for clusters)
            isInside = (ii_s>1) && (ii_e<S.Nx) && (jj_s>1) && (jj_e<S.Ny) && (kk_s>1) && (kk_e<S.Nz);
            assert(isInside,'Error: Atom too close to boundary for b calculation');
        elseif(S.BC == 2)
            % For periodic systems, find the overlap of rb-region and the fund. domain
            ii_s = max(ii_s,1);
            ii_e = min(ii_e,S.Nx);
            jj_s = max(jj_s,1);
            jj_e = min(jj_e,S.Ny);
            kk_s = max(kk_s,1);
            kk_e = min(kk_e,S.Nz);
        elseif(S.BC == 3)
            ii_s = max(ii_s,1);
            ii_e = min(ii_e,S.Nx);
            jj_s = max(jj_s,1);
            jj_e = min(jj_e,S.Ny);
            isInside = (kk_s>1) && (kk_e<S.Nz);
            assert(isInside,'Error: Atom too close to boundary in z-direction for b calculation');
        elseif(S.BC == 4)
            isInside = (ii_s>1) && (ii_e<S.Nx) && (jj_s>1) && (jj_e<S.Ny);
            assert(isInside,'Error: Atom too close to boundary for b calculation');
            kk_s = max(kk_s,1);
            kk_e = min(kk_e,S.Nz);
        end
        
        % % Check if the b-region is inside the domain (for cluster systems)
        % isInside = (ii_s>1) && (ii_e<S.Nx) && (jj_s>1) && (jj_e<S.Ny) && (kk_s>1) && (kk_e<S.Nz);
        % assert(isInside,'Error: Atom too close to boundary for b calculation');
        
        % Distances of grid points from atom JJ_a
        % count3 = 1;
        % dd = zeros(ii_e-ii_s+1+4*S.FDn, jj_e-jj_s+1+4*S.FDn, kk_e-kk_s+1+4*S.FDn);
        % for kk=kk_s-2*S.FDn:kk_e+2*S.FDn
        % zz = (kk - 1) * S.dz;
        % count2 = 1;
        % for jj=jj_s-2*S.FDn:jj_e+2*S.FDn
        % yy = (jj - 1) * S.dy;
        % count1 = 1;
        % for ii=ii_s-2*S.FDn:ii_e+2*S.FDn
        % xx = (ii - 1) * S.dx;
        % dd(count1,count2,count3) = norm([xx-x0, yy-y0, zz-z0]);
        % count1 = count1 + 1;
        % end
        % count2 = count2 + 1;
        % end
        % count3 = count3 + 1;
        % end
        % dd = reshape(dd,[],1);
        xx = (ii_s-2*S.FDn-1:ii_e+2*S.FDn-1)*S.dx - x0_i;
        yy = (jj_s-2*S.FDn-1:jj_e+2*S.FDn-1)*S.dy - y0_i;
        zz = (kk_s-2*S.FDn-1:kk_e+2*S.FDn-1)*S.dz - z0_i;
        [XX_3D,YY_3D,ZZ_3D] = meshgrid(xx,yy,zz);
        XX_3D = permute(XX_3D,[2 1 3]);
        YY_3D = permute(YY_3D,[2 1 3]);
        ZZ_3D = permute(ZZ_3D,[2 1 3]);
        % XX = reshape(XX_3D,[],1);
        % YY = reshape(YY_3D,[],1);
        % ZZ = reshape(ZZ_3D,[],1);
        dd = sqrt(S.metric_T(1,1)*XX_3D.^2 + S.metric_T(1,2)*(XX_3D.*YY_3D) + S.metric_T(1,3)*(XX_3D.*ZZ_3D) + ...
            S.metric_T(2,1)*(YY_3D.*XX_3D) + S.metric_T(2,2)*YY_3D.^2 + S.metric_T(2,3)*(YY_3D.*ZZ_3D) + ...
            S.metric_T(3,1)*(ZZ_3D.*XX_3D) + S.metric_T(3,2)*(ZZ_3D.*YY_3D) + S.metric_T(3,3)*ZZ_3D.^2) ;
        
        % Pseudopotential at grid points through interpolation
        %V_PS = interp1(r_grid_vloc, r_grid_vloc.*Vloc, dd, 'spline');
        
        %V_PS = interp1(r_grid_vloc, r_grid_vloc.*Vloc, dd, 'spline');
        V_PS = zeros(size(dd));
        IsLargeThanRmax = dd > r_grid_vloc(end);
        V_PS(IsLargeThanRmax) = -S.Atm(count_typ).Z;
        V_PS(~IsLargeThanRmax) = interp1(r_grid_vloc, r_grid_vloc.*Vloc, dd(~IsLargeThanRmax), 'spline');
        
        
        
        
        
        V_PS = V_PS./dd;
        V_PS(dd<r_grid_vloc(2)) = Vloc(1);
        % V_PS = reshape(V_PS,count1-1,count2-1,count3-1);
        % V_PS = reshape(V_PS,size(XX_3D));
        % Reference potential at grid points
        rc_ref = S.rc_ref; % WARNING: Might need smaller if pseudocharges overlap
        V_PS_ref = zeros(size(dd));
        I_ref = dd<rc_ref;
        V_PS_ref(~I_ref) = -(S.Atm(count_typ).Z)./dd(~I_ref);
        V_PS_ref(I_ref) = -S.Atm(count_typ).Z*(9*dd(I_ref).^7-30*rc_ref*dd(I_ref).^6 ...
            +28*rc_ref*rc_ref*dd(I_ref).^5-14*(rc_ref^5)*dd(I_ref).^2+12*rc_ref^7)/(5*rc_ref^8);
        % V_PS_ref = reshape(V_PS_ref,count1-1,count2-1,count3-1);
        % V_PS_ref = reshape(V_PS_ref,size(XX_3D));
        
        % Pseudocharge density
        % bJ = zeros(size(V_PS));
        % bJ_ref = zeros(size(V_PS));
        % count3 = S.FDn + 1;
        % for kk=kk_s-S.FDn:kk_e+S.FDn
        % count2 = S.FDn+1;
        % for jj=jj_s-S.FDn:jj_e+S.FDn
        % count1 = S.FDn+1;
        % for ii=ii_s-S.FDn:ii_e+S.FDn
        % bJ(count1,count2,count3) = (S.w2(1)/(dx2) + S.w2(1)/(dy2) + S.w2(1)/(dz2)) * V_PS(count1,count2,count3);
        % bJ_ref(count1,count2,count3) = (S.w2(1)/(dx2) + S.w2(1)/(dy2) + S.w2(1)/(dz2)) * V_PS_ref(count1,count2,count3);
        % for p=1:S.FDn
        % bJ(count1,count2,count3) = bJ(count1,count2,count3) + S.w2(p+1)/(dx2) * V_PS(count1+p,count2,count3) + ...
        % + S.w2(p+1)/(dx2) * V_PS(count1-p,count2,count3) + ...
        % + S.w2(p+1)/(dy2) * V_PS(count1,count2+p,count3) + ...
        % + S.w2(p+1)/(dy2) * V_PS(count1,count2-p,count3) + ...
        % + S.w2(p+1)/(dz2) * V_PS(count1,count2,count3+p) + ...
        % + S.w2(p+1)/(dz2) * V_PS(count1,count2,count3-p);
        % bJ_ref(count1,count2,count3) = bJ_ref(count1,count2,count3) + S.w2(p+1)/(dx2) * V_PS_ref(count1+p,count2,count3) + ...
        % + S.w2(p+1)/(dx2) * V_PS_ref(count1-p,count2,count3) + ...
        % + S.w2(p+1)/(dy2) * V_PS_ref(count1,count2+p,count3) + ...
        % + S.w2(p+1)/(dy2) * V_PS_ref(count1,count2-p,count3) + ...
        % + S.w2(p+1)/(dz2) * V_PS_ref(count1,count2,count3+p) + ...
        % + S.w2(p+1)/(dz2) * V_PS_ref(count1,count2,count3-p);
        % end
        % count1 = count1 + 1;
        % end
        % count2 = count2 + 1;
        % end
        % count3 = count3 + 1;
        % end
        
        % Pseudocharge density
        bJ = zeros(size(V_PS));
        bJ_ref = zeros(size(V_PS));
        II = 1+S.FDn : size(V_PS,1)-S.FDn;
        JJ = 1+S.FDn : size(V_PS,2)-S.FDn;
        KK = 1+S.FDn : size(V_PS,3)-S.FDn;
        bJ(II,JJ,KK) = coeff * V_PS(II,JJ,KK);
        bJ_ref(II,JJ,KK) = coeff * V_PS_ref(II,JJ,KK);
        for p = 1:S.FDn
            bJ(II,JJ,KK) = bJ(II,JJ,KK) + S.w2(p+1)*S.lapc_T(1,1)/dx2 * (V_PS(II+p,JJ,KK) + V_PS(II-p,JJ,KK)) + ...
                S.w2(p+1)*S.lapc_T(2,2)/dy2 * (V_PS(II,JJ+p,KK) + V_PS(II,JJ-p,KK)) + ...
                S.w2(p+1)*S.lapc_T(3,3)/dz2 * (V_PS(II,JJ,KK+p) + V_PS(II,JJ,KK-p));
            bJ_ref(II,JJ,KK) = bJ_ref(II,JJ,KK) + S.w2(p+1)*S.lapc_T(1,1)/dx2 * (V_PS_ref(II+p,JJ,KK) + V_PS_ref(II-p,JJ,KK)) + ...
                S.w2(p+1)*S.lapc_T(2,2)/dy2 * (V_PS_ref(II,JJ+p,KK) + V_PS_ref(II,JJ-p,KK)) + ...
                S.w2(p+1)*S.lapc_T(3,3)/dz2 * (V_PS_ref(II,JJ,KK+p) + V_PS_ref(II,JJ,KK-p));
            for q = 1:S.FDn
                bJ(II,JJ,KK) = bJ(II,JJ,KK) + S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,2) + S.lapc_T(2,1))/dxdy * ( V_PS(II+q,JJ+p,KK) - ...
                    V_PS(II-q,JJ+p,KK) - V_PS(II+q,JJ-p,KK) + V_PS(II-q,JJ-p,KK) ) + ...
                    S.w1(p+1)*S.w1(q+1)*(S.lapc_T(2,3) + S.lapc_T(3,2))/dydz * ( V_PS(II,JJ+q,KK+p) - ...
                    V_PS(II,JJ-q,KK+p) - V_PS(II,JJ+q,KK-p) + V_PS(II,JJ-q,KK-p) ) + ...
                    S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,3) + S.lapc_T(3,1))/dzdx * ( V_PS(II+q,JJ,KK+p) - ...
                    V_PS(II-q,JJ,KK+p) - V_PS(II+q,JJ,KK-p) + V_PS(II-q,JJ,KK-p) ) ;
                bJ_ref(II,JJ,KK) = bJ_ref(II,JJ,KK) + S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,2) + S.lapc_T(2,1))/dxdy * ( V_PS_ref(II+q,JJ+p,KK) - ...
                    V_PS_ref(II-q,JJ+p,KK) - V_PS_ref(II+q,JJ-p,KK) + V_PS_ref(II-q,JJ-p,KK) ) + ...
                    S.w1(p+1)*S.w1(q+1)*(S.lapc_T(2,3) + S.lapc_T(3,2))/dydz * ( V_PS_ref(II,JJ+q,KK+p) - ...
                    V_PS_ref(II,JJ-q,KK+p) - V_PS_ref(II,JJ+q,KK-p) + V_PS_ref(II,JJ-q,KK-p) ) + ...
                    S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,3) + S.lapc_T(3,1))/dzdx * ( V_PS_ref(II+q,JJ,KK+p) - ...
                    V_PS_ref(II-q,JJ,KK+p) - V_PS_ref(II+q,JJ,KK-p) + V_PS_ref(II-q,JJ,KK-p) ) ;
            end
        end
        bJ = (-1/(4*pi))*bJ;
        bJ_ref = (-1/(4*pi))*bJ_ref;
        
        % Calculate the gradient of pseudocharges
        dbJ_x = zeros(size(V_PS)); dbJ_y = zeros(size(V_PS)); dbJ_z = zeros(size(V_PS));
        dbJ_ref_x = zeros(size(V_PS)); dbJ_ref_y = zeros(size(V_PS)); dbJ_ref_z = zeros(size(V_PS));
        dVJ_x = zeros(size(V_PS)); dVJ_y = zeros(size(V_PS)); dVJ_z = zeros(size(V_PS));
        dVJ_ref_x = zeros(size(V_PS)); dVJ_ref_y = zeros(size(V_PS)); dVJ_ref_z = zeros(size(V_PS));
        II = 1+2*S.FDn : size(V_PS,1)-2*S.FDn;
        JJ = 1+2*S.FDn : size(V_PS,2)-2*S.FDn;
        KK = 1+2*S.FDn : size(V_PS,3)-2*S.FDn;
        for p = 1:S.FDn
            dbJ_x(II,JJ,KK) = dbJ_x(II,JJ,KK) + S.w1(p+1)/S.dx*(bJ(II+p,JJ,KK)-bJ(II-p,JJ,KK));
            dbJ_y(II,JJ,KK) = dbJ_y(II,JJ,KK) + S.w1(p+1)/S.dy*(bJ(II,JJ+p,KK)-bJ(II,JJ-p,KK));
            dbJ_z(II,JJ,KK) = dbJ_z(II,JJ,KK) + S.w1(p+1)/S.dz*(bJ(II,JJ,KK+p)-bJ(II,JJ,KK-p));
            dbJ_ref_x(II,JJ,KK) = dbJ_ref_x(II,JJ,KK) + S.w1(p+1)/S.dx*(bJ_ref(II+p,JJ,KK)-bJ_ref(II-p,JJ,KK));
            dbJ_ref_y(II,JJ,KK) = dbJ_ref_y(II,JJ,KK) + S.w1(p+1)/S.dy*(bJ_ref(II,JJ+p,KK)-bJ_ref(II,JJ-p,KK));
            dbJ_ref_z(II,JJ,KK) = dbJ_ref_z(II,JJ,KK) + S.w1(p+1)/S.dz*(bJ_ref(II,JJ,KK+p)-bJ_ref(II,JJ,KK-p));
            dVJ_x(II,JJ,KK) = dVJ_x(II,JJ,KK) + S.w1(p+1)/S.dx*(V_PS(II+p,JJ,KK)-V_PS(II-p,JJ,KK));
            dVJ_y(II,JJ,KK) = dVJ_y(II,JJ,KK) + S.w1(p+1)/S.dy*(V_PS(II,JJ+p,KK)-V_PS(II,JJ-p,KK));
            dVJ_z(II,JJ,KK) = dVJ_z(II,JJ,KK) + S.w1(p+1)/S.dz*(V_PS(II,JJ,KK+p)-V_PS(II,JJ,KK-p));
            dVJ_ref_x(II,JJ,KK) = dVJ_ref_x(II,JJ,KK) + S.w1(p+1)/S.dx*(V_PS_ref(II+p,JJ,KK)-V_PS_ref(II-p,JJ,KK));
            dVJ_ref_y(II,JJ,KK) = dVJ_ref_y(II,JJ,KK) + S.w1(p+1)/S.dy*(V_PS_ref(II,JJ+p,KK)-V_PS_ref(II,JJ-p,KK));
            dVJ_ref_z(II,JJ,KK) = dVJ_ref_z(II,JJ,KK) + S.w1(p+1)/S.dz*(V_PS_ref(II,JJ,KK+p)-V_PS_ref(II,JJ,KK-p));
        end
        % Calculate local force and correction force components
        [II_rb,JJ_rb,KK_rb] = meshgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
        II_rb = permute(II_rb,[2,1,3]);
        JJ_rb = permute(JJ_rb,[2,1,3]);
        KK_rb = permute(KK_rb,[2,1,3]);
        Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;
        %force_local(JJ_a,1) = force_local(JJ_a,1) + sum(sum(sum(S.W(Rowcount_rb).*dbJ_x(II,JJ,KK).*(S.phi(Rowcount_rb)))));
        %force_local(JJ_a,2) = force_local(JJ_a,2) + sum(sum(sum(S.W(Rowcount_rb).*dbJ_y(II,JJ,KK).*(S.phi(Rowcount_rb)))));
        %force_local(JJ_a,3) = force_local(JJ_a,3) + sum(sum(sum(S.W(Rowcount_rb).*dbJ_z(II,JJ,KK).*(S.phi(Rowcount_rb)))));
        
        force_local(JJ_a,1) = force_local(JJ_a,1) - sum(sum(sum(S.W(Rowcount_rb).*bJ(II,JJ,KK).*(Dphi_x(Rowcount_rb)))));
        force_local(JJ_a,2) = force_local(JJ_a,2) - sum(sum(sum(S.W(Rowcount_rb).*bJ(II,JJ,KK).*(Dphi_y(Rowcount_rb)))));
        force_local(JJ_a,3) = force_local(JJ_a,3) - sum(sum(sum(S.W(Rowcount_rb).*bJ(II,JJ,KK).*(Dphi_z(Rowcount_rb)))));
        
        %         force_corr(JJ_a,1) = force_corr(JJ_a,1) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( dbJ_ref_x(II,JJ,KK) .* (S.V_c(Rowcount_rb)-V_PS_ref(II,JJ,KK)) +...
        %             + dbJ_x(II,JJ,KK) .* (S.V_c(Rowcount_rb)+V_PS(II,JJ,KK)) +...
        %             + (dVJ_ref_x(II,JJ,KK)-dVJ_x(II,JJ,KK)).*(S.b_ref(Rowcount_rb)+S.b(Rowcount_rb)) +...
        %             + bJ(II,JJ,KK).*dVJ_x(II,JJ,KK) - bJ_ref(II,JJ,KK).*dVJ_ref_x(II,JJ,KK) ) )));
        %         force_corr(JJ_a,2) = force_corr(JJ_a,2) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( dbJ_ref_y(II,JJ,KK) .* (S.V_c(Rowcount_rb)-V_PS_ref(II,JJ,KK)) +...
        %             + dbJ_y(II,JJ,KK) .* (S.V_c(Rowcount_rb)+V_PS(II,JJ,KK)) +...
        %             + (dVJ_ref_y(II,JJ,KK)-dVJ_y(II,JJ,KK)).*(S.b_ref(Rowcount_rb)+S.b(Rowcount_rb)) +...
        %             + bJ(II,JJ,KK).*dVJ_y(II,JJ,KK) - bJ_ref(II,JJ,KK).*dVJ_ref_y(II,JJ,KK) ) )));
        %         force_corr(JJ_a,3) = force_corr(JJ_a,3) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( dbJ_ref_z(II,JJ,KK) .* (S.V_c(Rowcount_rb)-V_PS_ref(II,JJ,KK)) +...
        %             + dbJ_z(II,JJ,KK) .* (S.V_c(Rowcount_rb)+V_PS(II,JJ,KK)) +...
        %             + (dVJ_ref_z(II,JJ,KK)-dVJ_z(II,JJ,KK)).*(S.b_ref(Rowcount_rb)+S.b(Rowcount_rb)) +...
        %             + bJ(II,JJ,KK).*dVJ_z(II,JJ,KK) - bJ_ref(II,JJ,KK).*dVJ_ref_z(II,JJ,KK) ) )));
        %         force_corr(JJ_a,1) = force_corr(JJ_a,1) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_x(II,JJ,KK)-dVJ_x(II,JJ,KK)) +...
        %             + (dbJ_ref_x(II,JJ,KK) + dbJ_x(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
        %         force_corr(JJ_a,2) = force_corr(JJ_a,2) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_y(II,JJ,KK)-dVJ_y(II,JJ,KK)) +...
        %             + (dbJ_ref_y(II,JJ,KK) + dbJ_y(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
        %         force_corr(JJ_a,3) = force_corr(JJ_a,3) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_z(II,JJ,KK)-dVJ_z(II,JJ,KK)) +...
        %             + (dbJ_ref_z(II,JJ,KK) + dbJ_z(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
        
         
         force_corr(JJ_a,1) = force_corr(JJ_a,1) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( -(Dpseudo_x(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) +...
            + (dbJ_ref_x(II,JJ,KK) + dbJ_x(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
         force_corr(JJ_a,2) = force_corr(JJ_a,2) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( -(Dpseudo_y(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) +...
            + (dbJ_ref_y(II,JJ,KK) + dbJ_y(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
         force_corr(JJ_a,3) = force_corr(JJ_a,3) + 0.5*sum(sum(sum(S.W(Rowcount_rb).*( -(Dpseudo_z(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) +...
            + (dbJ_ref_z(II,JJ,KK) + dbJ_z(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
    end
    % Check if same type of atoms are over
    if count_typ_atms == S.Atm(count_typ).n_atm_typ
        count_typ_atms = 1;
        count_typ = count_typ + 1;
    else
        count_typ_atms = count_typ_atms + 1;
    end
    
    
end % end of loop over atoms


%**********************************************************************
%*                   Calculate nonlocal atomic force                  *
%**********************************************************************
% Type-I
for JJ_a = 1:S.n_atm % loop over all atoms
    if S.Atom(JJ_a).psptyp == 0
        force_kpt = zeros(S.tnkpt,3);
        for kpt = 1:S.tnkpt
            %for kpt = 1:S.tnkpt
            kpt_vec = S.kptgrid(kpt,:);
            % Calculate gradient of psi
            
            Dpsi_x = DiscreteGradient(S,kpt_vec(1),1)*S.psi(:,:,kpt);
            Dpsi_y = DiscreteGradient(S,kpt_vec(2),2)*S.psi(:,:,kpt);
            Dpsi_z = DiscreteGradient(S,kpt_vec(3),3)*S.psi(:,:,kpt);
            
            % Calculate nonlocal components of the force acting on atom JJ_a
            
            for l = 0:S.Atom(JJ_a).lmax
                if l == S.Atom(JJ_a).lloc
                    continue;
                end
                for m = -l:l
                    integral_1 = zeros(1,S.Nev);
                    integral_2_x = zeros(1,S.Nev);
                    integral_2_y = zeros(1,S.Nev);
                    integral_2_z = zeros(1,S.Nev);
                    for img = 1:S.Atom(JJ_a).n_image_rc
                        phase_fac = (exp(1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
                        integral_1 = integral_1 + (S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos,:,kpt)) * conj(phase_fac);
                        integral_2_x = integral_2_x + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                        integral_2_y = integral_2_y + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                        integral_2_z = integral_2_z + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                    end
                    tf_x = dot(S.occ(:,kpt),real(integral_1.*integral_2_x));
                    tf_y = dot(S.occ(:,kpt),real(integral_1.*integral_2_y));
                    tf_z = dot(S.occ(:,kpt),real(integral_1.*integral_2_z));
                    force_kpt(kpt,:) = force_kpt(kpt,:) + S.Atom(JJ_a).gamma_Jl(l+1) * [tf_x tf_y tf_z];
                end
            end
            
        end
        force_nloc(JJ_a,1) = -4*dot(S.wkpt,force_kpt(:,1));
        force_nloc(JJ_a,2) = -4*dot(S.wkpt,force_kpt(:,2));
        force_nloc(JJ_a,3) = -4*dot(S.wkpt,force_kpt(:,3));
        
    elseif S.Atom(JJ_a).psptyp == 1
        
        force_kpt = zeros(S.tnkpt,3);
        for kpt = 1:S.tnkpt
            %for kpt = 1:S.tnkpt
            kpt_vec = S.kptgrid(kpt,:);
            % Calculate gradient of psi
            
            Dpsi_x = DiscreteGradient(S,kpt_vec(1),1)*S.psi(:,:,kpt);
            Dpsi_y = DiscreteGradient(S,kpt_vec(2),2)*S.psi(:,:,kpt);
            Dpsi_z = DiscreteGradient(S,kpt_vec(3),3)*S.psi(:,:,kpt);
            
            % Calculate nonlocal components of the force acting on atom JJ_a
            
            for l = 0:S.Atom(JJ_a).lmax
                if l == S.Atom(JJ_a).lloc
                    continue;
                end
                for m = -l:l
                    integral_11 = zeros(1,S.Nev);
                    integral_12_x = zeros(1,S.Nev);
                    integral_12_y = zeros(1,S.Nev);
                    integral_12_z = zeros(1,S.Nev);
                    integral_21 = zeros(1,S.Nev);
                    integral_22_x = zeros(1,S.Nev);
                    integral_22_y = zeros(1,S.Nev);
                    integral_22_z = zeros(1,S.Nev);
                    for img = 1:S.Atom(JJ_a).n_image_rc
                        phase_fac = (exp(1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
                        integral_11 = integral_11 + (S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat1(:,m+l+1) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos,:,kpt)) * conj(phase_fac);
                        integral_12_x = integral_12_x + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat1(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                        integral_12_y = integral_12_y + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat1(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                        integral_12_z = integral_12_z + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat1(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                    end
                    for img = 1:S.Atom(JJ_a).n_image_rc
                        phase_fac = (exp(1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
                        integral_21 = integral_21 + (S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat2(:,m+l+1) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos,:,kpt)) * conj(phase_fac);
                        integral_22_x = integral_22_x + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat2(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                        integral_22_y = integral_22_y + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat2(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                        integral_22_z = integral_22_z + (conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat2(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * (Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * (phase_fac);
                    end
                    tf1_x = dot(S.occ(:,kpt),real(integral_11.*integral_12_x ));
                    tf1_y = dot(S.occ(:,kpt),real(integral_11.*integral_12_y ));
                    tf1_z = dot(S.occ(:,kpt),real(integral_11.*integral_12_z ));
                    tf2_x = dot(S.occ(:,kpt),real(integral_21.*integral_22_x ));
                    tf2_y = dot(S.occ(:,kpt),real(integral_21.*integral_22_y ));
                    tf2_z = dot(S.occ(:,kpt),real(integral_21.*integral_22_z ));
                    force_kpt(kpt,:) = force_kpt(kpt,:) + S.Atom(JJ_a).gamma_Jl1(l+1) * [tf1_x tf1_y tf1_z] + S.Atom(JJ_a).gamma_Jl2(l+1) * [tf2_x tf2_y tf2_z];
                end
            end
            
        end
        force_nloc(JJ_a,1) = -4*dot(S.wkpt,force_kpt(:,1));
        force_nloc(JJ_a,2) = -4*dot(S.wkpt,force_kpt(:,2));
        force_nloc(JJ_a,3) = -4*dot(S.wkpt,force_kpt(:,3));
         
    end
end % end of loop over atoms

%Type - II(not yet written for oncv psp)

% count_typ = 1;
% count_typ_atms = 1;
% for JJ_a = 1:S.n_atm % loop over all atoms
%     force_kpt = zeros(S.tnkpt,3);
%     for kpt = 1:S.tnkpt
%         kpt_vec = S.kptgrid(kpt,:);
%
%         % Calculate nonlocal components of the force acting on atom JJ_a
%
%         for l = 0:S.Atom(count_typ).lmax
%             if l == S.Atom(count_typ).lloc
%                 continue;
%             end
%             for m = -l:l
%                 integral_1 = zeros(1,S.Nev);
%                 integral_2_x = zeros(1,S.Nev);
%                 integral_2_y = zeros(1,S.Nev);
%                 integral_2_z = zeros(1,S.Nev);
%                 for img = 1:S.Atom(JJ_a).n_image_rc
%                     phase_fac = (exp(1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
%                     integral_1 = integral_1 + (S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)).' * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos,:,kpt)) * conj(phase_fac);
%                     xr =(S.Atom(JJ_a).rcImage(img).rc_pos_ii-1)*S.dx - S.Atom(JJ_a).rcImage(img).coordinates(1) ;
%                     yr =(S.Atom(JJ_a).rcImage(img).rc_pos_jj-1)*S.dy - S.Atom(JJ_a).rcImage(img).coordinates(2) ;
%                     zr =(S.Atom(JJ_a).rcImage(img).rc_pos_kk-1)*S.dz - S.Atom(JJ_a).rcImage(img).coordinates(3) ;
%                     dd_nl = sqrt(S.metric_T(1,1)*xr.^2 + S.metric_T(1,2)*(xr.*yr) + S.metric_T(1,3)*(xr.*zr) + ...
%                                  S.metric_T(2,1)*(yr.*xr) + S.metric_T(2,2)*yr.^2 + S.metric_T(2,3)*(yr.*zr) + ...
%                                  S.metric_T(3,1)*(zr.*xr) + S.metric_T(3,2)*(zr.*yr) + S.metric_T(3,3)*zr.^2);
%
%                     xr = xr(dd_nl > 1e-10);
%                     yr = yr(dd_nl > 1e-10);
%                     zr = zr(dd_nl > 1e-10);
%                     pos = S.Atom(JJ_a).rcImage(img).rc_pos(dd_nl > 1e-10);
%                     dd_nl = dd_nl(dd_nl > 1e-10);
%
%                     UdV_Jl = interp1(r_grid_vloc,UdV(:,l+1),dd_nl,'spline');
%                     [datasites,IA,IC] = unique(dd_nl);
%                     UdV_J = UdV_Jl(IA);
%                     Chi_spline = spline(datasites,UdV_J);
%                     DChi_spline = ppval(fnder(Chi_spline,1),datasites);
%                     Chi = UdV_J(IC);
%                     DChi = DChi_spline(IC);
%                     Ylm = SphericalHarmonics(S.lat_uvec,xr,yr,zr,l,m,'real');
%                     [DYlm_dx,DYlm_dy,DYlm_dz] = DYlm(S.lat_uvec,xr,yr,zr,l,m);
%                     integral_2_x = integral_2_x + (conj((DChi.*xr.*Ylm)./dd_nl + Chi.*DYlm_dx) .* S.W(pos)).' * (S.psi(pos,:,kpt)) * (phase_fac);
%                     integral_2_y = integral_2_y + (conj((DChi.*yr.*Ylm)./dd_nl + Chi.*DYlm_dy) .* S.W(pos)).' * (S.psi(pos,:,kpt)) * (phase_fac);
%                     integral_2_z = integral_2_z + (conj((DChi.*zr.*Ylm)./dd_nl + Chi.*DYlm_dz) .* S.W(pos)).' * (S.psi(pos,:,kpt)) * (phase_fac);
%                 end
%                 tf_x = dot(S.occ(:,kpt),real(integral_1.*integral_2_x));
%                 tf_y = dot(S.occ(:,kpt),real(integral_1.*integral_2_y));
%                 tf_z = dot(S.occ(:,kpt),real(integral_1.*integral_2_z));
%                 force_kpt(kpt,:) = force_kpt(kpt,:) + S.Atom(count_typ).gamma_Jl(l+1) * [tf_x tf_y tf_z];
%             end
%         end
%
%     end
%     force_nloc(JJ_a,1) = 4*dot(S.wkpt,force_kpt(:,1));
%     force_nloc(JJ_a,2) = 4*dot(S.wkpt,force_kpt(:,2));
%     force_nloc(JJ_a,3) = 4*dot(S.wkpt,force_kpt(:,3));
%
%     % Check if same type of atoms are over
%     if count_typ_atms == S.Atm(count_typ).n_atm_typ
%         count_typ_atms = 1;
%         count_typ = count_typ + 1;
%     else
%         count_typ_atms = count_typ_atms + 1;
%     end
%
% end % end of loop over atoms

force_loc = force_local + force_corr
force_nloc

force = force_local + force_corr + force_nloc;
force = force*S.grad_T;

end
