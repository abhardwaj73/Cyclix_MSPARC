function stress = evaluateStress_Cychel(S)
% @ brief    Function to calculate stress in cychel system (O(N^3))
% @ authors
%          Abhiraj Sharma <asharma424@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @2016-2019 (c) Georgia Institute of Technology.
% @ references
%             "On the calculation of the stress tensor in real-space Kohn-Sham
%              density functional theory (Sharma et. al. 2018)"
%===============================================================================
stress = 0.0;

% Stress contribution from Kinetic component
ks = 1;
for spin = 1:S.nspin
    for kpt = 1:S.tnkpt
        kpt_vec = S.kptgrid(kpt,:);

        Dpsi_z = blochGradient(S,kpt_vec,3)*S.psi(:,:,ks);
        Dcpsi_z = conj(Dpsi_z);
        
        stress = stress + real(-S.occfac * S.wkpt(kpt) * (transpose(S.W) * ... 
                      (Dcpsi_z.*Dpsi_z)) * S.occ(:,ks));
        
        ks = ks + 1;           
    end
end


% Stress contribution from exchange-correlation and energy terms from electrostatics
Drho_z = S.grad_3 * S.rho;
           
if S.nspin == 1
    stress = stress + S.Exc - sum(S.W .* S.Vxc .* S.rho)   + ...
           - sum(S.W' * (S.dvxcdgrho .* Drho_z .* Drho_z)) + ...
             0.5 * sum( S.W .* (S.b - S.rho) .* S.phi ) - S.Eself + S.E_corr ;
    
else
    stress = stress + S.Exc - sum(S.W'*(S.Vxc.*S.rho(:,2:3)))  + ...
           - sum(S.W' * (S.dvxcdgrho .* Drho_z .* Drho_z)) + ...
             0.5 * sum( S.W .* (S.b - S.rho(:,1)) .* S.phi ) - S.Eself + S.E_corr ;
end

% Stress contribution from remaining terms in electrostatics
Dphi_z = S.grad_3 * S.phi;
stress = stress + (1/4/pi) * sum( Dphi_z.*Dphi_z.*S.W);

%psdfilepath = sprintf('%s/PseudopotFiles',S.inputfile_path);
count_typ = 1;
count_typ_atms = 1;
for JJ_a = 1:S.n_atm % loop over all the atoms
    % Atom position
    x0 = S.Atoms(JJ_a,1);
    y0 = S.Atoms(JJ_a,2);
    z0 = S.Atoms(JJ_a,3);
    % Note the S.dx, S.dy, S.dz terms are to ensure the image rb-region overlap w/ fund. domain
    if S.BCx == 0
        n_image_xl = floor((S.Atoms(JJ_a,1) + S.Atm(count_typ).rb_x)/S.L1);
        n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+S.Atm(count_typ).rb_x-S.dx)/S.L1);
    else
        n_image_xl = 0;
        n_image_xr = 0;
    end
    
    if S.BCy == 0
        n_image_yl = floor((S.Atoms(JJ_a,2) + S.Atm(count_typ).rb_y)/S.L2);
        n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+S.Atm(count_typ).rb_y-S.dy)/S.L2);
    else
        n_image_yl = 0;
        n_image_yr = 0;
    end
    
    if S.BCz == 0
        n_image_zl = floor((S.Atoms(JJ_a,3) + S.Atm(count_typ).rb_z)/S.L3);
        n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+S.Atm(count_typ).rb_z-S.dz)/S.L3);
    else
        n_image_zl = 0;
        n_image_zr = 0;
    end
    
    % Total No. of images of atom JJ_a (including atom JJ_a)
    n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1) * (n_image_zl+n_image_zr+1);
    % Find the coordinates for all the images
    xx_img = [-n_image_xl : n_image_xr] * S.L1 + x0;
    yy_img = [-n_image_yl : n_image_yr] * S.L2 + y0;
    zz_img = [-n_image_zl : n_image_zr] * S.L3 + z0;
    [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = ndgrid(xx_img,yy_img,zz_img);

    % Loop over all image(s) of atom JJ_a (including atom JJ_a)
    for count_image = 1:n_image_total

        % Atom position of the image
        x0_i = XX_IMG_3D(count_image);
        y0_i = YY_IMG_3D(count_image);
        z0_i = ZZ_IMG_3D(count_image);

        % Indices of closest grid point to atom
        pos_ii = round((x0_i-S.xin) / S.dx) + 1;
        pos_jj = round((y0_i-S.yin) / S.dy) + 1;
        pos_kk = round((z0_i-S.zin) / S.dz) + 1;

        % Starting and ending indices of b-region
        ii_s = pos_ii - ceil(S.Atm(count_typ).rb_x/S.dx+0.5);
        ii_e = pos_ii + ceil(S.Atm(count_typ).rb_x/S.dx+0.5);
        jj_s = pos_jj - ceil(S.Atm(count_typ).rb_y/S.dy+0.5);
        jj_e = pos_jj + ceil(S.Atm(count_typ).rb_y/S.dy+0.5);
        kk_s = pos_kk - ceil(S.Atm(count_typ).rb_z/S.dz+0.5);
        kk_e = pos_kk + ceil(S.Atm(count_typ).rb_z/S.dz+0.5);

        % Check if the b-region is inside the domain in Dirichlet BC
        % direction
        isInside = (S.BCx == 0 || (S.BCx == 1 && (ii_s>1) && (ii_e<S.Nx))) && ...
           (S.BCy == 0 || (S.BCy == 1 && (jj_s>1) && (jj_e<S.Ny))) && ...
           (S.BCz == 0 || (S.BCz == 1 && (kk_s>1) && (kk_e<S.Nz)));
        assert(isInside,'Error: Atom too close to boundary for b calculation');
        ii_s = max(ii_s,1);
        ii_e = min(ii_e,S.Nx);
        jj_s = max(jj_s,1);
        jj_e = min(jj_e,S.Ny);
        kk_s = max(kk_s,1);
        kk_e = min(kk_e,S.Nz);

        xx = S.xin + (ii_s-2*S.FDn-1:ii_e+2*S.FDn-1)*S.dx;% - x0_i;
        yy = S.yin + (jj_s-2*S.FDn-1:jj_e+2*S.FDn-1)*S.dy;% - y0_i;
        zz = S.zin + (kk_s-2*S.FDn-1:kk_e+2*S.FDn-1)*S.dz;% - z0_i;
        [XX_3D,YY_3D,ZZ_3D] = ndgrid(xx,yy,zz);

        % Find distances
        dd = calculateDistance(XX_3D,YY_3D,ZZ_3D,x0_i,y0_i,z0_i,S);

        % Pseudopotential at grid points through interpolation
        V_PS = zeros(size(dd));
        IsLargeThanRmax = dd > S.Atm(count_typ).r_grid_vloc(end);
        V_PS(IsLargeThanRmax) = -S.Atm(count_typ).Z;
        V_PS(~IsLargeThanRmax) = interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).r_grid_vloc.*S.Atm(count_typ).Vloc, dd(~IsLargeThanRmax), 'spline');

        V_PS = V_PS./dd;
        V_PS(dd<S.Atm(count_typ).r_grid_vloc(2)) = S.Atm(count_typ).Vloc(1); % WARNING

        % Reference potential at grid points
        rc_ref = S.rc_ref; % WARNING: Might need smaller if pseudocharges overlap
        V_PS_ref = zeros(size(dd));
        I_ref = dd<rc_ref;
        V_PS_ref(~I_ref) = -(S.Atm(count_typ).Z)./dd(~I_ref);
        V_PS_ref(I_ref) = -S.Atm(count_typ).Z*(9*dd(I_ref).^7-30*rc_ref*dd(I_ref).^6 ...
            +28*rc_ref*rc_ref*dd(I_ref).^5-14*(rc_ref^5)*dd(I_ref).^2+12*rc_ref^7)/(5*rc_ref^8);

        % Pseudocharge density
        II = 1+S.FDn : size(V_PS,1)-S.FDn;
        JJ = 1+S.FDn : size(V_PS,2)-S.FDn;
        KK = 1+S.FDn : size(V_PS,3)-S.FDn;
        
        % Calculate bJ and bJ_ref
        bJ = pseudochargeDensity_atom(V_PS,II,JJ,KK,xx(1),S);
        bJ_ref = pseudochargeDensity_atom(V_PS_ref,II,JJ,KK,xx(1),S);

        bJ = (-1/(4*pi))*bJ;
        bJ_ref = (-1/(4*pi))*bJ_ref;

        % Calculate the gradient of pseudocharges
        II = 1+2*S.FDn : size(V_PS,1)-2*S.FDn;
        JJ = 1+2*S.FDn : size(V_PS,2)-2*S.FDn;
        KK = 1+2*S.FDn : size(V_PS,3)-2*S.FDn;
        [~, ~, dbJ_3] = dpseudopot(bJ,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
        [~, ~, dbJ_ref_3] = dpseudopot(bJ_ref,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
        [~, ~, dVJ_3] = dpseudopot(V_PS,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
        [~, ~, dVJ_ref_3] = dpseudopot(V_PS_ref,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);


        % Calculate local stress and correction stress components
        [II_rb,JJ_rb,KK_rb] = ndgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
        Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;

        [~,~,z1] = ndgrid((ii_s-1:ii_e-1)*S.dx - x0_i,(jj_s-1:jj_e-1)*S.dy - y0_i,(kk_s-1:kk_e-1)*S.dz - z0_i) ;
        
        stress = stress + sum(sum(sum( dVJ_3(II,JJ,KK) .* z1 .* ( - 0.5*bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) )));

        stress = stress + sum(sum(sum( dbJ_3(II,JJ,KK) .* z1 .* ( S.phi(Rowcount_rb) - 0.5*V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));

        stress = stress + 0.5 * sum(sum(sum( ( dbJ_3(II,JJ,KK) .* ( S.V_c(Rowcount_rb) + V_PS(II,JJ,KK) ) + ...
            dbJ_ref_3(II,JJ,KK) .* ( S.V_c(Rowcount_rb) - V_PS_ref(II,JJ,KK) ) + ...
            dVJ_ref_3(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ_ref(II,JJ,KK)) - ...
            dVJ_3(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ(II,JJ,KK)) ) .* z1 .* S.W(Rowcount_rb) )));   

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
%*                   Stress contribution from nonlocal          *
%**********************************************************************

for ks = 1:S.tnkpt*S.nspin
    if ks <= S.tnkpt
        kpt = ks;
    else
        kpt = ks - S.tnkpt;
    end

    if (kpt(1) == 0 && kpt(2) == 0 && kpt(3) == 0)
        fac = 1;
    else
        fac = 1i;
    end

    kpt_vec = S.kptgrid(kpt,:);
    
    TDpsi_3 = blochGradient(S,kpt_vec,3)*S.psi(:,:,ks);
    


    for JJ_a = 1:S.n_atm % loop over all atoms
        integral_1 = zeros(S.Atom(JJ_a).angnum,S.Nev);
        integral_2_zz = zeros(S.Atom(JJ_a).angnum,S.Nev);
        
        Chi_X_mult1 = zeros(S.Atom(JJ_a).angnum,S.Nev);
        
        for img = 1:S.Atom(JJ_a).n_image_rc
            phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
            Chi_X_mult1 = Chi_X_mult1 + transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chi_mat), S.W(S.Atom(JJ_a).rcImage(img).rc_pos))) * S.psi(S.Atom(JJ_a).rcImage(img).rc_pos,:,ks) * phase_fac ;
        end
        
        stress = stress - S.occfac * S.wkpt(kpt) * transpose(S.Atom(JJ_a).gamma_Jl) * (Chi_X_mult1.*conj(Chi_X_mult1)) * S.occ(:,ks) ;
        
        for img = 1:S.Atom(JJ_a).n_image_rc
            phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
            ChiW = transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chi_mat), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
            integral_1 = integral_1 + conj(ChiW) * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos,:,ks)) * conj(phase_fac);
            z_1 =(S.Atom(JJ_a).rcImage(img).rc_pos_kk-1)*S.dz - S.Atom(JJ_a).rcImage(img).coordinates(3) ;
            
            integral_2_zz = integral_2_zz + ChiW * ...
                ((TDpsi_3(S.Atom(JJ_a).rcImage(img).rc_pos,:)).*repmat(z_1,1,S.Nev)) * phase_fac ;
            
        end
        tf_zz = transpose(S.Atom(JJ_a).gamma_Jl) * real(integral_1.*integral_2_zz) * S.occ(:,ks);
        stress = stress - 2 * S.occfac * S.wkpt(kpt) * tf_zz;
        
    end % end of loop over atoms
end

cell_measure = S.L3;

stress = stress * (2*pi/S.L2)/ cell_measure;

