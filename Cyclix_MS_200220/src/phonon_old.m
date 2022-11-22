function [S] = phonon_old(S)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%  Calculate the eigenvalues and eigenvectors of the matrix of second 
	% derivative of the energy wrt atomic positions
	%  This will generate a square matrix
	%  Steps:
	%  1. Solve the Sternheimer equations for Dpsi_dR in a
	%     self-consistent mannerhive
	%  2. Evaluate the expression for the second derivative of energy
	%  
	% Disclaimer:
	% 1. No spin currently
	% 2. Check for poisson calculator
	% 3. Use dx*dy*dz instead of S.W
	% 4. Remove those states from Sternheimer for which occ is zero
	% 5. In the end make drho and dveff consistent and consistent with the dpsi
    % 6. Check symmetry of kpts and reduce the cost (k+q and -k+q)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	fid = fopen(S.DDBfname,'w') ;
	assert(fid~=-1,'Error: Cannot open .DDB file %s',S.DDBfname);

	fprintf(fid,'***************************************************************************\n');
	fprintf(fid,'                            Phononic properties                            \n');
	fprintf(fid,'***************************************************************************\n');

	fclose(fid);

    S.tnqpt = 1;
    N_c = round(2*pi/S.L2);
  	S.qptgrid = [0 0 0];
    S.qptgrid = S.qptgrid.*[1 1 2*pi/S.L3];
	
    SCF_tol = 5e-6;
	Poisson_tol = SCF_tol * 1e-3;
	S.linSolv_tol = 1e-8;

	S.Ns_nzocc = sum(S.occ > 1e-5); % is a function of kpt
	S.Ns_nzocc_max = max(S.Ns_nzocc);

	TOL = 1e-8;

    for qpt = 1:S.tnqpt

        if( abs(S.qptgrid(qpt,1)) < TOL && abs(S.qptgrid(qpt,2)) < TOL && abs(S.qptgrid(qpt,3)) < TOL )
            S.Kdeltaq0 = 1;
        else
            S.Kdeltaq0 = 0;
        end
      
        temp_epsilon = eps; % include the right boundary k-points instead of left
        kqptgrid = (S.kptgrid + S.qptgrid(qpt,:))./[1 1 2*pi/S.L3];
		kqptgrid_x = mod(kqptgrid(:,1) + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
		kqptgrid_y = mod(kqptgrid(:,2), N_c);
		kqptgrid_z = mod(kqptgrid(:,3) + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
        S.kqptgrid = [kqptgrid_x kqptgrid_y kqptgrid_z].*[1 1 2*pi/S.L3];
        
        
        mkqptgrid = (-S.kptgrid + S.qptgrid(qpt,:))./[1 1 2*pi/S.L3];
		mkqptgrid_x = mod(mkqptgrid(:,1) + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
		mkqptgrid_y = mod(mkqptgrid(:,2) + N_c, N_c);
		mkqptgrid_z = mod(mkqptgrid(:,3) + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
        S.mkqptgrid = [mkqptgrid_x mkqptgrid_y mkqptgrid_z].*[1 1 2*pi/S.L3];
        
     
        EigVal_kq = S.EigVal;
        occ_kq = S.occ;

        % mkqptgrid = S.mkqptgrid
        % aa
        
        
        if S.Kdeltaq0 == 0
            EigVal_mkq = S.EigVal;
            occ_mkq = S.occ;
            % S.psi_kq = zeros(S.N,S.Nev,S.tnkpt);
            % S.psi_mkq = zeros(S.N,S.Nev,S.tnkpt);
            for kpt = 1:S.tnkpt
                chkmatch1 = abs(S.kqptgrid(kpt,1) - S.kptgrid(:,1)) < TOL & abs(S.kqptgrid(kpt,2) - S.kptgrid(:,2)) < TOL & abs(S.kqptgrid(kpt,3) - S.kptgrid(:,3)) < TOL;
                chkmatch2 = abs(S.kqptgrid(kpt,1) + S.kptgrid(:,1)) < TOL & (abs(S.kqptgrid(kpt,2) + S.kptgrid(:,2) - N_c) < TOL | abs(S.kqptgrid(kpt,2) + S.kptgrid(:,2)) < TOL) & abs(S.kqptgrid(kpt,3) + S.kptgrid(:,3)) < TOL;
                if sum(chkmatch1) == 1
                    EigVal_kq(:,kpt) = S.EigVal(:,chkmatch1 == 1);
                    occ_kq(:,kpt) = S.occ(:,chkmatch1 == 1);
                    %S.psi_kq(:,:,kpt) = S.psi(:,:,chkmatch1 == 1);
                elseif sum(chkmatch2) == 1
                    EigVal_kq(:,kpt) = S.EigVal(:,chkmatch2 == 1);
                    occ_kq(:,kpt) = S.occ(:,chkmatch2 == 1);
                    %S.psi_kq(:,:,kpt) = conj(S.psi(:,:,chkmatch2 == 1));
                else
                    % Extract additional states and their eigenvalues from eigensolver
                    %[EigVal_kq(:,kpt),S.psi_kq(:,:,kpt)] = eigSolve_k(S,S.kqptgrid(kpt,:));
                end
     
                chkmatch1 = abs(S.mkqptgrid(kpt,1) - S.kptgrid(:,1)) < TOL & abs(S.mkqptgrid(kpt,2) - S.kptgrid(:,2)) < TOL & abs(S.mkqptgrid(kpt,3) - S.kptgrid(:,3)) < TOL;
                chkmatch2 = abs(S.mkqptgrid(kpt,1) + S.kptgrid(:,1)) < TOL & (abs(S.mkqptgrid(kpt,2) + S.kptgrid(:,2) - N_c) < TOL | abs(S.mkqptgrid(kpt,2) + S.kptgrid(:,2)) < TOL) & abs(S.mkqptgrid(kpt,3) + S.kptgrid(:,3)) < TOL;
                if sum(chkmatch1) == 1
                    EigVal_mkq(:,kpt) = S.EigVal(:,chkmatch1 == 1);
                    occ_mkq(:,kpt) = S.occ(:,chkmatch1 == 1);
                    %S.psi_mkq(:,:,kpt) = S.psi(:,:,chkmatch1 == 1);
                elseif sum(chkmatch2) == 1
                    EigVal_mkq(:,kpt) = S.EigVal(:,chkmatch2 == 1);
                    occ_mkq(:,kpt) = S.occ(:,chkmatch2 == 1);
                    %S.psi_mkq(:,:,kpt) = conj(S.psi(:,:,chkmatch2 == 1));
                else
                    % Extract additional states and their eigenvalues from eigensolver
                    %[EigVal_mkq(:,kpt),S.psi_mkq(:,:,kpt)] = eigSolve_k(S,S.mkqptgrid(kpt,:));
                end
            end
        end
        
        % Store coefficients for rhs of linear solve
        %===========================================
        S.coeff_Pv_psi_kq = zeros(S.Nev,S.Nev,S.tnkpt);
        %S.Beta_ph = zeros(S.Nev,S.Nev,S.tnkpt); % method-3
        S.Degn_index = zeros(S.Nev,S.Nev,S.tnkpt);

        for kpt = 1:S.tnkpt
            for i = 1:S.Nev
                for j = 1:S.Nev
                    temp = S.EigVal(i,kpt) - EigVal_kq(j,kpt);
                    temp_occ = S.occ(i,kpt) - occ_kq(j,kpt);
                    theta_mn = 1/(1+exp(S.bet*temp));% method-3
                    theta_nm = 1 - theta_mn;% method-3
                    if(j~=i)
                    	if(abs(temp) < 1e-4 && abs(temp_occ) < 1e-4)
                        	S.coeff_Pv_psi_kq(i,j,kpt) = 1 - (1-S.occ(i,kpt))*S.bet/2;
                        	S.Degn_index(i,j,kpt) = 1;
                            S.Degn_index(i,i,kpt) = 1;
                        	%S.Beta_ph(i,j,kpt) = S.occ(i,kpt)*theta_nm + S.occ(j,kpt)*theta_mn - S.occ(i,kpt)*(1-S.occ(i,kpt))*S.bet/2;% method-3           
                        else
                        	S.coeff_Pv_psi_kq(i,j,kpt) = 1/temp;

                       %  	if S.occ(i,kpt) > 1e-5
                       %  		S.coeff_Pv_psi_kq(i,j,kpt) = 1;
                     		% else
                       %  		S.coeff_Pv_psi_kq(i,j,kpt) = 1/temp;
                       %  	end
                        	%S.Beta_ph(i,j,kpt) = S.occ(i,kpt)*theta_nm + S.occ(j,kpt)*theta_mn + (S.occ(i,kpt)-S.occ(j,kpt))*theta_mn/temp;% method-3
                    	end
                    else
                    	if S.Kdeltaq0 == 1
                    		S.coeff_Pv_psi_kq(i,j,kpt) = 1;
                    	else
                    		if(abs(temp) < 1e-4 && abs(temp_occ) < 1e-4)
                    			S.coeff_Pv_psi_kq(i,j,kpt) = 1 - (1-S.occ(i,kpt))*S.bet/2;
                    		else
                    			S.coeff_Pv_psi_kq(i,j,kpt) = 1/temp;
                    			% if S.occ(i,kpt) > 1e-5
	                      %   		S.coeff_Pv_psi_kq(i,j,kpt) = 1;
	                     	% 	else
	                      %   		S.coeff_Pv_psi_kq(i,j,kpt) = 1/temp;
	                      %   	end
                    		end
                    	end
                    end

                end
            end
            S.coeff_Pv_psi_kq(:,:,kpt) = transpose(S.coeff_Pv_psi_kq(:,:,kpt));
            %S.Beta_ph(:,:,kpt) = transpose(S.Beta_ph(:,:,kpt));% method-3
        end
        
        %a = S.coeff_Pv_psi_kq
        %aa
        % b = S.Degn_index
        
        
        if S.Kdeltaq0 == 0
            S.coeff_Pv_psi_mkq = zeros(S.Nev,S.Nev,S.tnkpt);
            for kpt = 1:S.tnkpt
                for i = 1:S.Nev
                    for j = 1:S.Nev
                        temp = S.EigVal(i,kpt) - EigVal_mkq(j,kpt);
                        temp_occ = S.occ(i,kpt) - occ_mkq(j,kpt);
                        theta_mn = 1/(1+exp(S.bet*temp));% method-3
                        theta_nm = 1 - theta_mn;% method-3
                        if(abs(temp) < 1e-4 && abs(temp_occ) < 1e-4)
                            S.coeff_Pv_psi_mkq(i,j,kpt) = 1 - (1-S.occ(i,kpt))*S.bet/2;
                        else
                        	S.coeff_Pv_psi_mkq(i,j,kpt) = 1/temp;
                       %      if S.occ(i,kpt) > 1e-5
                       %  		S.coeff_Pv_psi_mkq(i,j,kpt) = 1;
                     		% else
                       %  		S.coeff_Pv_psi_mkq(i,j,kpt) = 1/temp;
                       %  	end
                        end
                    end
                end
                S.coeff_Pv_psi_mkq(:,:,kpt) = transpose(S.coeff_Pv_psi_mkq(:,:,kpt));
            end
        end

       % a = S.coeff_Pv_psi_mkq

		% 1. Guess drho/dR
		%==================
		
		S.Drho_dR = zeros(S.N,3*S.n_atm);

		% Evaluate dphi/dR
		%==================
		S.M_dyn = zeros(3*S.n_atm,3*S.n_atm);
		S.M_dyn_var = zeros(3*S.n_atm,3*S.n_atm);
		S.M_dyn_xc = zeros(3*S.n_atm,3*S.n_atm);
		S.M_dyn_elec = zeros(3*S.n_atm,3*S.n_atm);
		S.M_dyn_nl = zeros(3*S.n_atm,3*S.n_atm);
		S.M_dyn_occ = zeros(3*S.n_atm,3*S.n_atm);
		S.Dphi_dR = zeros(S.N,3*S.n_atm);
		S.Db_dR = zeros(S.N,3*S.n_atm);
		%S = const_for_FFT(S);
		S = electrostatics_local(S,qpt);
		%S.Drho_dR = repmat(S.rho - sum(S.rho)/S.N,[1,3*S.n_atm]);

		% Evaluate dVxc_dR
		%==================
		S.Dvxc_dR = bsxfun(@times,S.Dvxc_drho,S.Drho_dR);

		% Evaluate dVeff_dR
		%===================
		S.Dveff_dR = bsxfun(@plus,S.Dphi_dR,S.Dvxc_dR);
	

		% SCF procedure for solving the Sternheimer equations
		% (Total of 3*n_atm scf solves)
		%=====================================================

		% No. of states to be considered for DFPT
		S.Dpsi_dR_kq = ones(S.N,S.Ns_nzocc_max,S.tnkpt);
        if S.Kdeltaq0 == 0
            S.Dpsi_dR_mkq = ones(S.N,S.Ns_nzocc_max,S.tnkpt);
        end
		S.Dlambda_dR = zeros(S.Ns_nzocc_max,S.tnkpt,3*S.n_atm);
		S.DlambdaF_dR = zeros(3*S.n_atm,1);
		S.Dg_dR = zeros(S.Ns_nzocc_max,S.tnkpt,3*S.n_atm);

		for scf_natmc = 1:3*S.n_atm
			% Initialize the mixing history vectors
			S.X = 0*S.X;
			S.F = 0*S.F;
			S.mixing_hist_fkm1 = 0*S.mixing_hist_fkm1;

			if S.MixingVariable == 0
				S.mixing_hist_xkm1 = S.Drho_dR(:,scf_natmc);
				rho_temp = S.Drho_dR(:,scf_natmc);
			else
				% for potential mixing, we store the mean-0 part only
				if S.BC == 2
					Veff_mean = mean(S.Dveff_dR(:,scf_natmc));
				else
					Veff_mean = 0.0;
				end
				S.mixing_hist_xkm1 = S.Dveff_dR(:,scf_natmc) - Veff_mean;
				Veff_temp = S.Dveff_dR(:,scf_natmc);
			end

			count_SCF = 1;
			err = 100;
			max_scf_iter = 100;
			min_scf_iter = 3;
			if max_scf_iter < min_scf_iter
				min_scf_iter = max_scf_iter;
			end

			% start scf loop
			while (err > SCF_tol && count_SCF <= max_scf_iter || count_SCF <= min_scf_iter)
				tic_scf = tic;
				fprintf(' ======================================================= \n');
				fprintf(' Sternheimer equation: %2d \n SCF iteration number: %2d \n',scf_natmc,count_SCF);
				fprintf(' ======================================================= \n');

				if S.parallel ~= 1
					S = sternheimer(S,scf_natmc,qpt);
				else
					S = sternheimer_kparal(S,scf_natmc,qpt);
				end

				% Solve for Fermi energy S.lambda_f and occupations
				if S.Kdeltaq0 == 1
					S = occupation_pert(S,scf_natmc); % method-3 off it
				end

				% Electron density
				S = electronDensity_pert(S,scf_natmc);

				if S.MixingVariable == 1
					% update Veff
					
					% Electrostatic potential
					rhs = S.Drho_dR(:,scf_natmc) + S.Db_dR(:,scf_natmc);
					f = poisson_RHS(S,rhs);
					[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,S.qptgrid(qpt,:));
					%Hfun = @(x) aar_phonon(S,[],[],[],[],DL11,DL22,DL33,DG1,DG2,DG3,[],f,x,[],Poisson_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,2);
					%S.Dphi_dR(:,scf_natmc) = gmres(Hfun,f,[],1e-6,1000);
					S.Dphi_dR(:,scf_natmc) = aar_phonon(S,[],[],[],[],DL11,DL22,DL33,DG1,DG2,DG3,[],f,S.Dphi_dR(:,scf_natmc),[],Poisson_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,2);
					%S.Dphi_dR(:,scf_natmc) = poissonSolve_FFT(S,rhs);

					if(S.BC == 2 && S.Kdeltaq0 == 1)
						% To make sure integral phi is 0 (this removes the arbitariness of the
						% constant in phi calculation)
						S.Dphi_dR(:,scf_natmc) = S.Dphi_dR(:,scf_natmc) - dot(S.W,S.Dphi_dR(:,scf_natmc))/sum(S.W); % change weight here
                    end

					% Exchange-correlation potential
					S.Dvxc_dR(:,scf_natmc) = bsxfun(@times,S.Dvxc_drho,S.Drho_dR(:,scf_natmc));

					% Effective potential
					S.Dveff_dR(:,scf_natmc) = bsxfun(@plus,S.Dphi_dR(:,scf_natmc),S.Dvxc_dR(:,scf_natmc));
		        end
				
				% Error in SCF fixed-point iteration
				if S.MixingVariable == 1
					err = (norm(S.Dveff_dR(:,scf_natmc) - Veff_temp))/(norm(S.Dveff_dR(:,scf_natmc)));
				else
					err = (norm(S.Drho_dR(:,scf_natmc) - rho_temp))/(norm(S.Drho_dR(:,scf_natmc)));
				end
					
				fprintf(' Error in SCF iteration: %.4e \n',err) ;

				% Mixing to accelerate SCF convergence
				if (err > SCF_tol || count_SCF < S.MINIT_SCF)
					if S.MixingVariable == 1 % potential mixing
						% shift Veff and Veff_temp so they have mean 0
						if S.BC == 2 && S.Kdeltaq0 == 1
							Veff_mean = mean(S.Dveff_dR(:,scf_natmc)); Veff_temp_mean = mean(Veff_temp);
						else
							Veff_mean = 0; Veff_temp_mean = 0;
						end
						S.Dveff_dR(:,scf_natmc) = S.Dveff_dR(:,scf_natmc) - Veff_mean;
						Veff_temp = Veff_temp - Veff_temp_mean;
						% enter mixing
						[S,S.Dveff_dR(:,scf_natmc)] = mixing(S,S.Dveff_dR(:,scf_natmc),Veff_temp,count_SCF);
						% shift the mean back
						S.Dveff_dR(:,scf_natmc) = S.Dveff_dR(:,scf_natmc) + Veff_mean; % the new veff for next input
						% note we add Veff_mean not Veff_temp_mean, since the one
						% used for the next input is S.Veff 
						Veff_temp = S.Dveff_dR(:,scf_natmc);
					else % density mixing
						[S, S.Drho_dR(:,scf_natmc)] = mixing(S,S.Drho_dR(:,scf_natmc),rho_temp,count_SCF);
						rho_temp = S.Drho_dR(:,scf_natmc);
						% at this point rho_temp = S.rho, i.e., new input density
						% update Veff
						rhs = S.Drho_dR(:,scf_natmc) + S.Db_dR(:,scf_natmc);
						f = poisson_RHS(S,rhs);
						[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,S.qptgrid(qpt,:));
						%Hfun = @(x) aar_phonon(S,[],[],[],[],DL11,DL22,DL33,DG1,DG2,DG3,[],f,x,[],Poisson_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,2);
						%S.Dphi_dR(:,scf_natmc) = gmres(Hfun,f,[],1e-6,1000);
						S.Dphi_dR(:,scf_natmc) = aar_phonon(S,[],[],[],[],DL11,DL22,DL33,DG1,DG2,DG3,[],f,S.Dphi_dR(:,scf_natmc),[],Poisson_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,2);
					
						%S.Dphi_dR(:,scf_natmc) = aar(S.Lap_std,f,S.Dphi_dR(:,scf_natmc),Poisson_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU);
						%S.Dphi_dR(:,scf_natmc) = poissonSolve_FFT(S,rhs);

						if(S.BC == 2 && S.Kdeltaq0 == 1 )
							% To make sure integral phi is 0 (this removes the arbitariness of the
							% constant in phi calculation)
							S.Dphi_dR(:,scf_natmc) = S.Dphi_dR(:,scf_natmc) - dot(S.W,S.Dphi_dR(:,scf_natmc))/sum(S.W);
						end

						% Exchange-correlation potential
						S.Dvxc_dR(:,scf_natmc) = bsxfun(@times,S.Dvxc_drho,S.Drho_dR(:,scf_natmc));

						% Effective potential
						S.Dveff_dR(:,scf_natmc) = bsxfun(@plus,S.Dphi_dR(:,scf_natmc),S.Dvxc_dR(:,scf_natmc));
					end
						
				end
				count_SCF = count_SCF + 1;
				scf_runtime = toc(tic_scf);
				fprintf(' This SCF iteration took %.3f s.\n\n', scf_runtime);
			end

			if (count_SCF == max_scf_iter + 1)
				disp(' SCF did not converge. Maximum iterations reached!')
			end

			fprintf('\n Finished SCF iteration in %d steps!\n', (count_SCF - 1));
			fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

% 			if scf_natmc < 3*S.n_atm
% 				S.Drho_dR(:,scf_natmc+1) = S.Drho_dR(:,scf_natmc);
% 				S.Dveff_dR(:,scf_natmc+1) = S.Dveff_dR(:,scf_natmc);
% 			end

			S = Calculate_dynamical_matrix(S,scf_natmc);

		end


      	if S.Kdeltaq0 == 1
      		M_dyn_var_lt = 2 * real(S.M_dyn_var);
        	S.M_dyn_var = M_dyn_var_lt + transpose(M_dyn_var_lt) - diag(diag(M_dyn_var_lt));
        	S.M_dyn_occ = S.M_dyn_occ + transpose(S.M_dyn_occ) - diag(diag(S.M_dyn_occ));
        else
        	S.M_dyn_var = S.M_dyn_var + transpose(conj(S.M_dyn_var)) - diag(diag(S.M_dyn_var));
        end

		S.M_dyn_xc = transpose(conj(S.Drho_dR)) * bsxfun(@times,S.Dvxc_dR, S.W);
		S.M_dyn_elec = S.M_dyn_elec + 0.5*(transpose(conj(S.Drho_dR+S.Db_dR)) * bsxfun(@times,S.Dphi_dR, S.W) + transpose(conj(S.Dphi_dR)) * bsxfun(@times,(S.Drho_dR+S.Db_dR), S.W)); 
		wkpt = S.wkpt;

		for kpt = 1:S.tnkpt
			kpt_vec = S.kptgrid(kpt,:);
			fac = 1.0i;
			
			PsiW = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
			
			for JJ_a = 1:S.n_atm % loop over all atoms

				% Calculate nonlocal components of the force acting on atom JJ_a
				Chi_full_k = zeros(S.N,S.Atom(JJ_a).angnum);
				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
		 			phase_fac_k = exp(-1*dot(kpt_vec,img_disp*fac));
		 			Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_k;
				end

				DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
		        DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
		        DChi_z_temp_k = blochGradient(S,S.kptgrid(kpt,:),3) * Chi_full_k;
		        DChi_x = zeros(S.N,S.Atom(JJ_a).angnum);
		        DChi_y = zeros(S.N,S.Atom(JJ_a).angnum);
		        for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
				DChi_z = DChi_z_temp_k;
		        
				DChi_xx_temp_k = blochGradient(S,kpt_vec,1)*DChi_x;
				DChi_yx_temp_k = blochGradient(S,kpt_vec,2)*DChi_x;
				DChi_zx_temp_k = blochGradient(S,kpt_vec,3)*DChi_x;
				DChi_xy_temp_k = blochGradient(S,kpt_vec,1)*DChi_y;
				DChi_yy_temp_k = blochGradient(S,kpt_vec,2)*DChi_y;
				DChi_zy_temp_k = blochGradient(S,kpt_vec,3)*DChi_y;
				DChi_xz_temp_k = blochGradient(S,kpt_vec,1)*DChi_z;
				DChi_yz_temp_k = blochGradient(S,kpt_vec,2)*DChi_z;				
				DChi_zz_temp_k = blochGradient(S,kpt_vec,3)*DChi_z;
				DChi_xx = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_yx = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_zx = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_xy = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_yy = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_zy = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_xz = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_yz = zeros(S.N,S.Atom(JJ_a).angnum);
				DChi_zz = DChi_zz_temp_k;

				% integral_xx = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_yx = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_zx = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_xy = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_yy = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_zy = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_xz = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_yz = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);
				% integral_zz = zeros(S.Ns_nzocc(kpt),S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					
					DChi_xx(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_xx(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_xx_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_yx_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_xy(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_xy(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_xy_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_yy_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_xz(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_xz(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_xz_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_yz_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);

					DChi_yx(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_yx(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_xx_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_yx_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_yy(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_yy(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_xy_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_yy_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_yz(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_yz(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_xz_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_yz_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);																																												
					
					DChi_zx(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_zx(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(3,3)*DChi_zx_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_zy(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_zy(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(3,3)*DChi_zy_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				

					%ROTi = ROT';
					% integral_xx = integral_xx - ROT(1,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,1) - ROT(1,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,1) - ...
					%                             ROT(1,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,1) - ROT(1,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,1);

					% integral_xy = integral_xy - ROT(1,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,2) - ROT(1,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,2) - ...
					%                             ROT(1,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,2) - ROT(1,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,2);

					% integral_xz = integral_xz - ROT(1,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_z_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(3,3) - ROT(1,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_z_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(3,3);

					% integral_yx = integral_yx - ROT(2,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,1) - ROT(2,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,1) - ...
					%                             ROT(2,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,1) - ROT(2,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,1);

					% integral_yy = integral_yy - ROT(2,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,2) - ROT(2,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,2) - ...
					%                             ROT(2,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,2) - ROT(2,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,2);                           

					% integral_yz = integral_yz - ROT(2,1) * dpsiW_x(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_z_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(3,3) - ROT(2,2) * dpsiW_y(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_z_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(3,3);

					% integral_zx = integral_zx - ROT(3,3) * dpsiW_z(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,1) - ROT(3,3) * dpsiW_z(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,1);

					% integral_zy = integral_zy - ROT(3,3) * dpsiW_z(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(1,2) - ROT(3,3) * dpsiW_z(:,S.Atom(JJ_a).rcImage(img).rc_pos) * DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) * ROTi(2,2);



				end

				%integral_zz = -dpsiW_z *  DChi_z_temp_k;


				% integral_xx
				% integral_yx
				% integral_zx
				% integral_xy
				% integral_yy
				% integral_zy

				% integral_xz
				% integral_yz

				

				
				integral = conj(PsiW) * conj(Chi_full_k);
				integral_x = -1*PsiW*DChi_x;
				integral_y = -1*PsiW*DChi_y;
				integral_z = -1*PsiW*DChi_z;

				% integral_xT = -1*PsiW*DChi_xT;
				% integral_yT = -1*PsiW*DChi_yT;
				
				integral_xx = PsiW*DChi_xx;
				integral_yx = PsiW*DChi_yx;
				integral_zx = PsiW*DChi_zx;
				integral_xy = PsiW*DChi_xy;
				integral_yy = PsiW*DChi_yy;
				integral_zy = PsiW*DChi_zy;

				integral_xz = PsiW*DChi_xz;
				integral_yz = PsiW*DChi_yz;
				integral_zz = PsiW*DChi_zz;

				S.M_dyn_nl(3*(JJ_a-1)+1,3*(JJ_a-1)+1) = S.M_dyn_nl(3*(JJ_a-1)+1,3*(JJ_a-1)+1) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_x.*conj(integral_x)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_xx.*integral) * S.Atom(JJ_a).gamma_Jl);
		        S.M_dyn_nl(3*(JJ_a-1)+1,3*(JJ_a-1)+2) = S.M_dyn_nl(3*(JJ_a-1)+1,3*(JJ_a-1)+2) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_x.*conj(integral_y)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_xy.*integral) * S.Atom(JJ_a).gamma_Jl);														  
				S.M_dyn_nl(3*(JJ_a-1)+1,3*(JJ_a-1)+3) = S.M_dyn_nl(3*(JJ_a-1)+1,3*(JJ_a-1)+3) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_x.*conj(integral_z)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_xz.*integral) * S.Atom(JJ_a).gamma_Jl);																			
				S.M_dyn_nl(3*(JJ_a-1)+2,3*(JJ_a-1)+1) = S.M_dyn_nl(3*(JJ_a-1)+2,3*(JJ_a-1)+1) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_y.*conj(integral_x)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_yx.*integral) * S.Atom(JJ_a).gamma_Jl);
		        S.M_dyn_nl(3*(JJ_a-1)+2,3*(JJ_a-1)+2) = S.M_dyn_nl(3*(JJ_a-1)+2,3*(JJ_a-1)+2) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_y.*conj(integral_y)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_yy.*integral) * S.Atom(JJ_a).gamma_Jl);														  
				S.M_dyn_nl(3*(JJ_a-1)+2,3*(JJ_a-1)+3) = S.M_dyn_nl(3*(JJ_a-1)+2,3*(JJ_a-1)+3) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_y.*conj(integral_z)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_yz.*integral) * S.Atom(JJ_a).gamma_Jl);																			
				S.M_dyn_nl(3*(JJ_a-1)+3,3*(JJ_a-1)+1) = S.M_dyn_nl(3*(JJ_a-1)+3,3*(JJ_a-1)+1) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_z.*conj(integral_x)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_zx.*integral) * S.Atom(JJ_a).gamma_Jl);
		        S.M_dyn_nl(3*(JJ_a-1)+3,3*(JJ_a-1)+2) = S.M_dyn_nl(3*(JJ_a-1)+3,3*(JJ_a-1)+2) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_z.*conj(integral_y)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_zy.*integral) * S.Atom(JJ_a).gamma_Jl);														  
				S.M_dyn_nl(3*(JJ_a-1)+3,3*(JJ_a-1)+3) = S.M_dyn_nl(3*(JJ_a-1)+3,3*(JJ_a-1)+3) + S.occfac * wkpt(kpt) * ...
																							(transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_z.*conj(integral_z)) * S.Atom(JJ_a).gamma_Jl + ...
																							 transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * real(integral_zz.*integral) * S.Atom(JJ_a).gamma_Jl);

			end
	    end

	    if S.Kdeltaq0 == 1
        	S.M_dyn_nl = real(S.M_dyn_nl) + transpose(real(S.M_dyn_nl));
        else
        	S.M_dyn_nl = S.M_dyn_nl + transpose(conj(S.M_dyn_nl));
        end

		% Full dynamical matrix
		S.M_dyn = S.M_dyn_var + S.M_dyn_xc + S.M_dyn_elec + S.M_dyn_nl + S.M_dyn_occ;
		

		fileID = fopen(S.DDBfname,'a');
		if (fileID == -1) 
			error('\n Cannot open file "%s"\n',S.DDBfname);
		end
		fprintf(fileID, '\nDynamical matrix before symmetrization (not shifted) (Ha/Bohr^2) for q = [%18.10f %18.10f %18.10f]:\n',S.qptgrid(1),S.qptgrid(2),S.qptgrid(3));
		for i = 1:3*S.n_atm
			for j = 1:3*S.n_atm
				fprintf(fileID, '%18.10f + %18.10fi\t', real(S.M_dyn(i,j)), imag(S.M_dyn(i,j)));
			end
			fprintf(fileID,'\n');
		end
		fclose(fileID);

		S.M_dyn = (S.M_dyn + S.M_dyn')/2;


		fileID = fopen(S.DDBfname,'a');
		if (fileID == -1) 
			error('\n Cannot open file "%s"\n',S.DDBfname);
		end
		fprintf(fileID, '\nDynamical matrix after symmetrization (not shifted) (Ha/Bohr^2) for q = [%18.10f %18.10f %18.10f]:\n',S.qptgrid(1),S.qptgrid(2),S.qptgrid(3));
		for i = 1:3*S.n_atm
			for j = 1:3*S.n_atm
				fprintf(fileID, '%18.10f + %18.10fi\t', real(S.M_dyn(i,j)), imag(S.M_dyn(i,j)));
			end
			fprintf(fileID,'\n');
		end
		fclose(fileID);

		%S.M_dyn = (S.M_dyn + transpose(S.M_dyn))/2;

		% Acoustic sum rule (ASR)
		%=========================
		if S.Kdeltaq0 == 1 && S.BC <= 2
			for row = 1:S.n_atm
				for indx = 1:3
					S.M_dyn(3*(row-1)+indx,3*(row-1)+1) = -sum(S.M_dyn(3*(row-1)+indx,1:3:3*(S.n_atm-1)+1)) + S.M_dyn(3*(row-1)+indx,3*(row-1)+1);
					S.M_dyn(3*(row-1)+indx,3*(row-1)+2) = -sum(S.M_dyn(3*(row-1)+indx,2:3:3*(S.n_atm-1)+2)) + S.M_dyn(3*(row-1)+indx,3*(row-1)+2);
					S.M_dyn(3*(row-1)+indx,3*(row-1)+3) = -sum(S.M_dyn(3*(row-1)+indx,3:3:3*(S.n_atm-1)+3)) + S.M_dyn(3*(row-1)+indx,3*(row-1)+3);
				end
			end


			fileID = fopen(S.DDBfname,'a');
			if (fileID == -1) 
				error('\n Cannot open file "%s"\n',S.DDBfname);
			end
			fprintf(fileID, '\nDynamical matrix after shift (Ha/Bohr^2) for q = [%18.10f %18.10f %18.10f]:\n',S.qptgrid(1),S.qptgrid(2),S.qptgrid(3));
			for i = 1:3*S.n_atm
				for j = 1:3*S.n_atm
					fprintf(fileID, '%18.10f + %18.10fi\t', real(S.M_dyn(i,j)), imag(S.M_dyn(i,j)));
				end
				fprintf(fileID,'\n');
			end
			fclose(fileID);
		end


		% Phonon modes and frequencies
		%=============================
		M_mass = S.M_mass * S.amu2au;
		S.M_IFC = M_mass\S.M_dyn;

		[S.phonon_mode, S.phonon_omega2] = eig(S.M_IFC);

		S.phonon_omega = sqrt(S.phonon_omega2) * S.fs2atu * 1e5/(2.99792458)/2/pi; % in cm-1

		fileID = fopen(S.DDBfname,'a');
		if (fileID == -1) 
			error('\n Cannot open file "%s"\n',S.DDBfname);
		end
		fprintf(fileID, '\nPhonon_frequencies (cm-1) for q = [%18.10f %18.10f %18.10f]:\n',S.qptgrid(1),S.qptgrid(2),S.qptgrid(3));
		for i = 1:3*S.n_atm
			fprintf(fileID, '%18.10f + %18.10fi\t', real(S.phonon_omega(i,i)), imag(S.phonon_omega(i,i))) ;
		end

		fprintf(fileID, '\n Phonon modes for q = [%18.10f %18.10f %18.10f]:\n',S.qptgrid(1),S.qptgrid(2),S.qptgrid(3));
		for i = 1:3*S.n_atm
			for j = 1:3*S.n_atm
				fprintf(fileID, '%18.10f\t', S.phonon_mode(i,j));
			end
			fprintf(fileID,'\n');
		end
		fclose(fileID);

		%M_dyn = S.M_dyn

		%phonon_omega = S.phonon_omega
    end

end






function [S] = Calculate_dynamical_matrix(S,scf_natmc)

	wkpt = S.wkpt;
        

	Dlambda_dR = zeros(S.Ns_nzocc_max,3*S.n_atm);

	for atmc = 1:scf_natmc
		Dir = rem(atmc,3);
		if Dir == 0
			Dir = 3;
		end

		for kpt = 1:S.tnkpt
			kqpt_vec = S.kqptgrid(kpt,:);				
			fac = 1.0i;
            
            % Method-2
			%==============================================================================================================

			JJ_a = ceil(atmc/3);
			Chi_full_k = zeros(S.N,S.Atom(JJ_a).angnum);
			Chi_full_kq = zeros(S.N,S.Atom(JJ_a).angnum);
			for img = 1:S.Atom(JJ_a).n_image_rc
				img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	 			phase_fac_k = exp(-1*dot(S.kptgrid(kpt,:),img_disp*fac));
	 			phase_fac_kq = exp(-1*dot(kqpt_vec,img_disp*fac));
				Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_k;
				Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_kq;
			end

			if Dir == 3
	        	DChi_k = blochGradient(S,S.kptgrid(kpt,:),Dir) * Chi_full_k;
	        	DChi_kq = blochGradient(S,S.kqptgrid(kpt,:),Dir) * Chi_full_kq;
	        elseif Dir == 1
	        	DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
	        	DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
	        	DChi_x_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),1) * Chi_full_kq;
	        	DChi_y_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),2) * Chi_full_kq;
	         	DChi_k = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_kq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			else
				DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
	        	DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
	        	DChi_x_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),1) * Chi_full_kq;
	        	DChi_y_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),2) * Chi_full_kq;
	         	DChi_k = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_kq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			end
			
		    atmc2 = scf_natmc;
	    	PsinW = transpose(bsxfun(@times, S.psi(:,1:S.Ns_nzocc(kpt),kpt), S.W));
			DVeff_psin = bsxfun(@times,S.Dveff_dR(:,atmc),S.psi(:,1:S.Ns_nzocc(kpt),kpt));

		    integral_1 = PsinW * conj(Chi_full_k);
		    integral_2 = -1*PsinW * conj(DChi_k);

		    DH_psin = DVeff_psin + Chi_full_kq * bsxfun(@times,transpose(integral_2),S.Atom(JJ_a).gamma_Jl) - ...
				      DChi_kq * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);

			Dlambda_dR(1:S.Ns_nzocc(kpt),atmc2) = diag(conj(PsinW) * DH_psin);
            S.M_dyn_var(atmc2,atmc) = S.M_dyn_var(atmc2,atmc) - 2 * wkpt(kpt) * S.occ(1:S.Ns_nzocc(kpt),kpt)' * (diag(bsxfun(@times,S.Dpsi_dR_kq(:,1:S.Ns_nzocc(kpt),kpt),S.W)' * DH_psin));

			%=======================================================================================================================

			temp = S.occ(1:S.Ns_nzocc(kpt),kpt) .* (1-S.occ(1:S.Ns_nzocc(kpt),kpt));
			ns_nzocc = (temp > 1e-4);
			if (sum (ns_nzocc) > 0)
				atmc2 = scf_natmc;
                Dg_dR = (-S.bet) * temp(ns_nzocc) .*(Dlambda_dR(ns_nzocc,atmc2) - S.DlambdaF_dR(atmc));
				S.M_dyn_occ(atmc2,atmc) = S.M_dyn_occ(atmc2,atmc) + 2 * wkpt(kpt) * (1/S.bet) * sum(S.Dg_dR(ns_nzocc,kpt,atmc2) .* Dg_dR ./temp(ns_nzocc));% - ...
																    %2 * S.wkpt(kpt) * S.DlambdaF_dR(atmc2) * sum(S.Dg_dR(ns_nzocc,kpt,atmc)) - 2 * S.wkpt(kpt) * S.DlambdaF_dR(atmc) * sum(S.Dg_dR(ns_nzocc,kpt,atmc2));
			end

		end
    end
        
        
    if S.Kdeltaq0 == 0
        for atmc = 1:scf_natmc
            Dir = rem(atmc,3);
            if Dir == 0
                Dir = 3;
            end
            
            for kpt = 1:S.tnkpt    
                mkqpt_vec = S.mkqptgrid(kpt,:);
                fac = 1.0i;

                % Method-2
                %==============================================================================================================
                
                JJ_a = ceil(atmc/3);
                
                Chi_full_mk = zeros(S.N,S.Atom(JJ_a).angnum);
                Chi_full_mkq = zeros(S.N,S.Atom(JJ_a).angnum);
                for img = 1:S.Atom(JJ_a).n_image_rc
                    img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
                    phase_fac_mk = exp(-1*dot(-S.kptgrid(kpt,:),img_disp*fac));
                    phase_fac_mkq = exp(-1*dot(mkqpt_vec,img_disp*fac));
                    Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mk;
                    Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mkq;
                end
                
               	if Dir == 3
		        	DChi_mk = blochGradient(S,-S.kptgrid(kpt,:),Dir) * Chi_full_mk;
		        	DChi_mkq = blochGradient(S,S.mkqptgrid(kpt,:),Dir) * Chi_full_mkq;
		        elseif Dir == 1
		        	DChi_x_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),1) * Chi_full_mk;
		        	DChi_y_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),2) * Chi_full_mk;
		        	DChi_x_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),1) * Chi_full_mkq;
		        	DChi_y_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),2) * Chi_full_mkq;
		         	DChi_mk = zeros(S.N,S.Atom(JJ_a).angnum);
		         	DChi_mkq = zeros(S.N,S.Atom(JJ_a).angnum);

					for img = 1:S.Atom(JJ_a).n_image_rc
						img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
						fac1 = -img_disp(2)/S.L2;
						fac2 = -img_disp(3)/S.L3;
						ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
						DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
						DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					end
				else
					DChi_x_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),1) * Chi_full_mk;
		        	DChi_y_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),2) * Chi_full_mk;
		        	DChi_x_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),1) * Chi_full_mkq;
		        	DChi_y_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),2) * Chi_full_mkq;
		         	DChi_mk = zeros(S.N,S.Atom(JJ_a).angnum);
		         	DChi_mkq = zeros(S.N,S.Atom(JJ_a).angnum);

					for img = 1:S.Atom(JJ_a).n_image_rc
						img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
						fac1 = -img_disp(2)/S.L2;
						fac2 = -img_disp(3)/S.L3;
						ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
						DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
						DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					end
				end
                
                
                atmc2 = scf_natmc;
                PsiW_mk = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
                DVeff_psi = bsxfun(@times,S.Dveff_dR(:,atmc),conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)));
                
                integral_1 = PsiW_mk * conj(Chi_full_mk);
                integral_2 = -1*PsiW_mk * conj(DChi_mk);
                
                DH_psi = DVeff_psi + Chi_full_mkq * bsxfun(@times,transpose(integral_2),S.Atom(JJ_a).gamma_Jl) - ...
                         DChi_mkq * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
                
                S.M_dyn_var(atmc2,atmc) = S.M_dyn_var(atmc2,atmc) - 2 * wkpt(kpt) * S.occ(1:S.Ns_nzocc(kpt),kpt)' * (diag(bsxfun(@times,S.Dpsi_dR_mkq(:,1:S.Ns_nzocc(kpt),kpt),S.W)' * DH_psi));
            end
        end
    end
            
        
        

	

	for kpt = 1:S.tnkpt
		kpt_vec = S.kptgrid(kpt,:);
		kqpt_vec = S.kqptgrid(kpt,:);
		fac = 1.0i;
		
		PsiW = transpose(bsxfun(@times, (S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
		DPsiW = transpose(bsxfun(@times, conj(S.Dpsi_dR_kq(:,1:S.Ns_nzocc(kpt),kpt)), S.W));

		for JJ_a = 1:S.n_atm % loop over all atoms

			% Calculate nonlocal components of the force acting on atom JJ_a
			Chi_full_k = zeros(S.N,S.Atom(JJ_a).angnum);
			Chi_full_kq = zeros(S.N,S.Atom(JJ_a).angnum);
			for img = 1:S.Atom(JJ_a).n_image_rc
				img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	 			phase_fac_k = exp(-1*dot(kpt_vec,img_disp*fac));
	 			phase_fac_kq = exp(-1*dot(kqpt_vec,img_disp*fac));
				Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_k;
				Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_kq;
			end

			DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
	        DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
	        DChi_z_temp_k = blochGradient(S,S.kptgrid(kpt,:),3) * Chi_full_k;
	        DChi_x_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),1) * Chi_full_kq;
	        DChi_y_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),2) * Chi_full_kq;
	        DChi_z_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),3) * Chi_full_kq;
	        DChi_x = zeros(S.N,S.Atom(JJ_a).angnum);
	        DChi_kqx = zeros(S.N,S.Atom(JJ_a).angnum);
	        DChi_y = zeros(S.N,S.Atom(JJ_a).angnum);
	        DChi_kqy = zeros(S.N,S.Atom(JJ_a).angnum);
	        for img = 1:S.Atom(JJ_a).n_image_rc
				img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
				fac1 = -img_disp(2)/S.L2;
				fac2 = -img_disp(3)/S.L3;
				ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
				DChi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				DChi_kqx(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kqx(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				DChi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				DChi_kqy(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kqy(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
			end
			DChi_z = DChi_z_temp_k;
	        DChi_kqz = DChi_z_temp_kq;
	        
	        Psi_Chistar = PsiW*conj(Chi_full_k);

			Psi_Dchistar_x = -1*PsiW*conj(DChi_x);
			Psi_Dchistar_y = -1*PsiW*conj(DChi_y);
			Psi_Dchistar_z = -1*PsiW*conj(DChi_z);

				

			Dpsistar_Chi = DPsiW*Chi_full_kq;
			Dpsistar_Dchi_x = -1*DPsiW*DChi_kqx;
			Dpsistar_Dchi_y = -1*DPsiW*DChi_kqy;
			Dpsistar_Dchi_z = -1*DPsiW*DChi_kqz;
			


			S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+1) = S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+1) + (S.Kdeltaq0+1) * S.occfac * wkpt(kpt) * ...
										         (S.Kdeltaq0 * transpose(S.Dg_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc)) * (conj(Psi_Dchistar_x).*Psi_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
										          transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Dchi_x.*Psi_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
										   		  transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Chi.*Psi_Dchistar_x) * S.Atom(JJ_a).gamma_Jl);
			S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+2) = S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+2) + (S.Kdeltaq0+1) * S.occfac * wkpt(kpt) * ...
										         (S.Kdeltaq0 * transpose(S.Dg_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc)) * (conj(Psi_Dchistar_y).*Psi_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
										          transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Dchi_y.*Psi_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
										   		  transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Chi.*Psi_Dchistar_y) * S.Atom(JJ_a).gamma_Jl);
			S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+3) = S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+3) + (S.Kdeltaq0+1) * S.occfac * wkpt(kpt) * ...
										         (S.Kdeltaq0 * transpose(S.Dg_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc)) * (conj(Psi_Dchistar_z).*Psi_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
										          transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Dchi_z.*Psi_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
										   		  transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Chi.*Psi_Dchistar_z) * S.Atom(JJ_a).gamma_Jl);							         
		end
    end
        
    if S.Kdeltaq0 == 0
        for kpt = 1:S.tnkpt
            mkpt_vec = -S.kptgrid(kpt,:);
            mkqpt_vec = S.mkqptgrid(kpt,:);
            fac = 1.0i;
            
            PsiW = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
		    DPsiW = transpose(bsxfun(@times, conj(S.Dpsi_dR_mkq(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
            
            for JJ_a = 1:S.n_atm % loop over all atoms
                
                % Calculate nonlocal components of the force acting on atom JJ_a
                Chi_full_mk = zeros(S.N,S.Atom(JJ_a).angnum);
                Chi_full_mkq = zeros(S.N,S.Atom(JJ_a).angnum);
                for img = 1:S.Atom(JJ_a).n_image_rc
                    img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
                    phase_fac_mk = exp(-1*dot(mkpt_vec,img_disp*fac));
                    phase_fac_mkq = exp(-1*dot(mkqpt_vec,img_disp*fac));
                    Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mk;
                    Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mkq;
                end
            

                DChi_x_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),1) * Chi_full_mk;
		        DChi_y_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),2) * Chi_full_mk;
		        DChi_z_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),3) * Chi_full_mk;
		        DChi_x_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),1) * Chi_full_mkq;
		        DChi_y_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),2) * Chi_full_mkq;
		        DChi_z_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),3) * Chi_full_mkq;
		        DChi_x = zeros(S.N,S.Atom(JJ_a).angnum);
		        DChi_kqx = zeros(S.N,S.Atom(JJ_a).angnum);
		        DChi_y = zeros(S.N,S.Atom(JJ_a).angnum);
		        DChi_kqy = zeros(S.N,S.Atom(JJ_a).angnum);
		        for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_x(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kqx(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kqx(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_y(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kqy(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kqy(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
				DChi_z = DChi_z_temp_mk;
		        DChi_kqz = DChi_z_temp_mkq;
		        
		        Psistar_Chistar = PsiW*conj(Chi_full_mk);

				Psistar_Dchistar_x = -1*PsiW*conj(DChi_x);
				Psistar_Dchistar_y = -1*PsiW*conj(DChi_y);
				Psistar_Dchistar_z = -1*PsiW*conj(DChi_z);			

				Dpsistar_Chi = DPsiW*Chi_full_mkq;
				Dpsistar_Dchi_x = -1*DPsiW*DChi_kqx;
				Dpsistar_Dchi_y = -1*DPsiW*DChi_kqy;
				Dpsistar_Dchi_z = -1*DPsiW*DChi_kqz;  
				
                
                S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+1) = S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+1) + S.occfac * wkpt(kpt) * ...
							                         (transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Dchi_x.*Psistar_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
							                         transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Chi.*Psistar_Dchistar_x) * S.Atom(JJ_a).gamma_Jl);
				S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+2) = S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+2) + S.occfac * wkpt(kpt) * ...
							                         (transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Dchi_y.*Psistar_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
							                         transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Chi.*Psistar_Dchistar_y) * S.Atom(JJ_a).gamma_Jl);
							                         			                         
                S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+3) = S.M_dyn_nl(scf_natmc,3*(JJ_a-1)+3) + S.occfac * wkpt(kpt) * ...
							                         (transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Dchi_z.*Psistar_Chistar) * S.Atom(JJ_a).gamma_Jl + ...
							                         transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)) * (Dpsistar_Chi.*Psistar_Dchistar_z) * S.Atom(JJ_a).gamma_Jl);
        	end
    	end
    end
    
        
        
end












function S = electronDensity_pert(S,scf_natmc)

	S.Drho_dR(:,scf_natmc) = 0;
	if S.Kdeltaq0 == 1
		for kpt =1:S.tnkpt
			S.Drho_dR(:,scf_natmc) = S.Drho_dR(:,scf_natmc) + 2 * S.wkpt(kpt) * (sum( bsxfun(@times,S.psi(:,1:S.Ns_nzocc(kpt),kpt).*conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)),S.Dg_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc)'),2) + ...
		 													  2 * real(sum( bsxfun(@times,S.Dpsi_dR_kq(:,1:S.Ns_nzocc(kpt),kpt).*conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)),S.occ(1:S.Ns_nzocc(kpt),kpt)'),2))) ;
			%S.Drho_dR(:,scf_natmc) = S.Drho_dR(:,scf_natmc) + 2 * S.wkpt(kpt) * (2 * real(sum(S.Dpsi_dR_kq(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc).*conj(S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc)),2))) ; method-3
		end
	else
		for kpt =1:S.tnkpt
			S.Drho_dR(:,scf_natmc) = S.Drho_dR(:,scf_natmc) + 2 * S.wkpt(kpt) * (sum( bsxfun(@times,S.Dpsi_dR_mkq(:,1:S.Ns_nzocc(kpt),kpt).*S.psi(:,1:S.Ns_nzocc(kpt),kpt),S.occ(1:S.Ns_nzocc(kpt),kpt)'),2) + ...
		 													                     sum( bsxfun(@times,S.Dpsi_dR_kq(:,1:S.Ns_nzocc(kpt),kpt).*conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)),S.occ(1:S.Ns_nzocc(kpt),kpt)'),2));
		end
	end
	sum(S.Drho_dR(:,scf_natmc).*S.W)
	%dot(S.Drho_dR(:,scf_natmc),S.Drho_dR(:,scf_natmc).*S.W)

end



function S = occupation_pert(S,scf_natmc)

	num = 0;
	denom = 0;
	for kpt = 1:S.tnkpt
		G_times_one_minus_G = S.occ(1:S.Ns_nzocc(kpt),kpt) .* (1-S.occ(1:S.Ns_nzocc(kpt),kpt)); 
		num = num + S.wkpt(kpt) * dot(G_times_one_minus_G,S.Dlambda_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc));
		denom = denom + S.wkpt(kpt) * sum(G_times_one_minus_G);
	end

	if denom < 1e-4
		S.DlambdaF_dR(scf_natmc) = 0;
		S.Dg_dR(:,:,scf_natmc) = 0;
	else
		S.DlambdaF_dR(scf_natmc) = num/denom;
		S.Dg_dR(:,:,scf_natmc) = (-S.bet) * S.occ(1:S.Ns_nzocc_max,:) .* (1-S.occ(1:S.Ns_nzocc_max,:)) .*(S.Dlambda_dR(:,:,scf_natmc) - S.DlambdaF_dR(scf_natmc)) ;
	end

	%S.Dg_dR(:,:,scf_natmc)

end




function S = sternheimer(S,scf_natmc,qpt)

	rhs_strnhimr = zeros(S.N,S.Ns_nzocc_max);

	Dir = rem(scf_natmc,3);
	if Dir == 0
	    Dir = 3;
	end

	if S.Kdeltaq0 == 1
		for kpt = 1:S.tnkpt
	        %--------------------------Method-1--------------------------------------------------
	        %====================================================================================
	        % PsiW = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
	        % DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi(:,1:S.Ns_nzocc(kpt),kpt));
	        % psi_DVeff_psi = PsiW * DVeff_psi;
	        
	        % % non-local psp (can be done only once during a scf)
	        % kpt_vec = S.kptgrid(kpt,:);
	        % fac = 1.0i;
	        % JJ_a = ceil(scf_natmc/3);
	        % Chi_full = zeros(S.N,S.Atom(JJ_a).angnum);
	        % for img = 1:S.Atom(JJ_a).n_image_rc
	        % 	img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	        %	phase_fac = exp(-1*dot(kpt_vec,img_disp*fac));
	        % 	Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac;
	        % end
	        % DChi = blochGradient(S,kpt_vec,Dir) * Chi_full;
	        
	        % chi_star_psi = conj(PsiW) * conj(Chi_full);
	        % psi_star_Dchi = -1*PsiW*DChi;
	        
	        % temp =  psi_star_Dchi * bsxfun(@times,transpose(chi_star_psi),S.Atom(JJ_a).gamma_Jl);
	        % psi_DVnl_psi = temp + temp';
	        
	        
	        
	        % % For degenerate case: STARTS
	        % %=============================
	        % psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        % %psi_DH_psi = 0.5 * (psi_DH_psi + psi_DH_psi')
	        % M_Pertb_degn =  psi_DH_psi .* S.Degn_index(1:S.Ns_nzocc(kpt),1:S.Ns_nzocc(kpt),kpt); %0.5 * (psi_DH_psi + psi_DH_psi')
	        % [CC,CW] = eig(M_Pertb_degn,'nobalance');
	        
	        % S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt) * CC;
	        
	        % PsinW = transpose(bsxfun(@times, conj(S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc)), S.W));
	        % DVeff_psin = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc));
	        % psin_DVeff_psin = PsinW * DVeff_psin;
	        % psi_DVeff_psin = PsiW * DVeff_psin;
	        
	        % chi_star_psin = conj(PsinW) * conj(Chi_full);
	        % psin_star_Dchi = -1*PsinW*DChi;
	        % temp =  psin_star_Dchi * bsxfun(@times,transpose(chi_star_psin),S.Atom(JJ_a).gamma_Jl);
	        % psin_DVnl_psin = temp + temp';
	        
	        % psi_star_chi = conj(chi_star_psi);
	        % Dchi_star_psin = conj(psin_star_Dchi);
	        
	        % temp1 = psi_star_chi * bsxfun(@times,transpose(Dchi_star_psin),S.Atom(JJ_a).gamma_Jl);
	        % temp2 = psi_star_Dchi * bsxfun(@times,transpose(chi_star_psin),S.Atom(JJ_a).gamma_Jl);
	        % psi_DVnl_psin = temp1 + temp2;
	        
	        % % For degenerate case: ENDS
	        % %============================
	        
	        % psin_DH_psin = psin_DVeff_psin + psin_DVnl_psin;
	        % psi_DH_psin = psi_DVeff_psin + psi_DVnl_psin;
	        % S.Dlambda_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc) = diag(psin_DH_psin);
	        % rhs_strnhimr(:,1:S.Ns_nzocc(kpt),kpt) = bsxfun(@times,S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc), transpose(S.Dlambda_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc))) + ...
	        %                                         S.psi(:,1:S.Ns_nzocc(kpt),kpt) * (psi_DH_psin .* S.coeff_Pv_psi(1:S.Ns_nzocc(kpt),1:S.Ns_nzocc(kpt),kpt)) - ...
	        % 				                        Chi_full * bsxfun(@times,transpose(Dchi_star_psin),S.Atom(JJ_a).gamma_Jl) + ...
	        % 									    DChi * bsxfun(@times,transpose(chi_star_psin),S.Atom(JJ_a).gamma_Jl) - DVeff_psin;
	        
	        %-----------------------------------Method-1 ends--------------------------------------------
	        
	        %-----------------------------------Method-2-------------------------------------------------
	        %============================================================================================

	        PsiW = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
	        DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi(:,1:S.Ns_nzocc(kpt),kpt));
	        psi_DVeff_psi = PsiW * DVeff_psi;

	        % non-local psp (can be done only once during a scf)
	        kpt_vec = S.kptgrid(kpt,:);
	        fac = 1.0i;
	        JJ_a = ceil(scf_natmc/3);
	        Chi_full = zeros(S.N,S.Atom(JJ_a).angnum);
	        for img = 1:S.Atom(JJ_a).n_image_rc
	            img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	            phase_fac = exp(-1*dot(kpt_vec,img_disp*fac));
	            Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac;
	        end

	        if Dir == 3
	        	DChi = blochGradient(S,kpt_vec,Dir) * Chi_full;
	        elseif Dir == 1
	        	DChi_x_temp = blochGradient(S,kpt_vec,1)*Chi_full;
	 			DChi_y_temp = blochGradient(S,kpt_vec,2)*Chi_full;
	         	DChi = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			else
				DChi_x_temp = blochGradient(S,kpt_vec,1)*Chi_full;
	 			DChi_y_temp = blochGradient(S,kpt_vec,2)*Chi_full;
	         	DChi = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			end





	        
	        integral_1 = conj(PsiW) * conj(Chi_full);
	        integral_2 = -1*PsiW*DChi;
	        
	        temp =  integral_2 * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
	        psi_DVnl_psi = temp + temp';
	        
	        %For degenerate case: STARTS
	        %=============================
	        % psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        % M_Pertb_degn = psi_DH_psi .* S.Degn_index(1:S.Ns_nzocc(kpt),1:S.Ns_nzocc(kpt),kpt);
	        % %M_Pertb_degn = 0.5*(M_Pertb_degn + M_Pertb_degn');
	        % [CC,~] = eig(M_Pertb_degn,'nobalance');
	        % %[CC(2:4,2:4),DD] = eig(M_Pertb_degn(2:4,2:4));
	        % %[CC(5:7,5:7),DD] = eig(M_Pertb_degn(5:7,5:7));
	        
	        % S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt) * (CC);
	        
	        % PsiW = transpose(bsxfun(@times, conj(S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc)), S.W));
	        % DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc));
	        % psi_DVeff_psi = PsiW * DVeff_psi;
	        
	        % integral_1 = conj(PsiW) * conj(Chi_full);
	        % integral_2 = -1*PsiW*DChi;
	        % temp =  integral_2 * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
	        % psi_DVnl_psi = temp + temp';
	        %For degenerate case: ENDS
	        %============================

	        %S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt);
	        
	        psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        S.Dlambda_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc) = diag(psi_DH_psi);
	        coeff_Pv_Psi = S.coeff_Pv_psi_kq(1:S.Ns_nzocc(kpt),1:S.Ns_nzocc(kpt),kpt);% + diag(ones(S.Ns_nzocc(kpt),1));
	        rhs_strnhimr(:,1:S.Ns_nzocc(kpt)) = S.psi(:,1:S.Ns_nzocc(kpt),kpt) * (psi_DH_psi .* coeff_Pv_Psi) - Chi_full * bsxfun(@times,integral_2',S.Atom(JJ_a).gamma_Jl) + ...
	                                            DChi * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl) - DVeff_psi;

	  %      	rhs_strnhimr(26385:26395,1,kpt)
			% rhs_strnhimr(52235:52245,1,kpt)
			% aa       
	        %---------------------------Method-2 ends-----------------------------------------------------
	        
	        %---------------------------Method-3 (Barroni et. al.)----------------------------------------
	        %=============================================================================================
    	% 	PsiW = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
    	% 	DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi(:,1:S.Ns_nzocc(kpt),kpt));
    	% 	psi_DVeff_psi = PsiW * DVeff_psi;
    
    	% 	% non-local psp (can be done only once during a scf)
    	% 	kpt_vec = S.kptgrid(kpt,:);
    	% 	fac = 1.0i;
    	% 	JJ_a = ceil(scf_natmc/3);
    	% 	Chi_full = zeros(S.N,S.Atom(JJ_a).angnum);
    	% 	for img = 1:S.Atom(JJ_a).n_image_rc
    	% 		img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
     % 			phase_fac = exp(-1*dot(kpt_vec,img_disp*fac));
    	% 		Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac;
    	% 	end
    	% 	DChi = blochGradient(S,kpt_vec,Dir) * Chi_full;
    
    	% 	integral_1 = conj(PsiW) * conj(Chi_full);
    	% 	integral_2 = -1*PsiW*DChi;
    
    	% 	temp =  integral_2 * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
    	% 	psi_DVnl_psi = temp + temp';
    
    	% 	% %For degenerate case: STARTS
    	% 	% %=============================
    	% 	% psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
    	% 	% M_Pertb_degn = psi_DH_psi .* S.Degn_index(1:S.Ns_nzocc(kpt),1:S.Ns_nzocc(kpt),kpt);
    	% 	% [CC,~] = eig(M_Pertb_degn,'nobalance');
    	% 	% S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt) * CC;
    
    	% 	% PsiW = transpose(bsxfun(@times, conj(S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc)), S.W));
    	% 	% DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc));
    	% 	% psi_DVeff_psi = PsiW * DVeff_psi;
    
    	% 	% integral_1 = conj(PsiW) * conj(Chi_full);
    	% 	% integral_2 = -1*PsiW*DChi;
    	% 	% temp =  integral_2 * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
    	% 	% psi_DVnl_psi = temp + temp';
    	% 	% %For degenerate case: ENDS
    	% 	% %============================
    	% 	S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt);		
    
    	% 	psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
    	% 	S.Dlambda_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc) = diag(psi_DH_psi);
    
    	% 	rhs_strnhimr(:,1:S.Ns_nzocc(kpt),kpt) = S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) * (psi_DH_psi .* S.Beta_ph(1:S.Ns_nzocc(kpt),1:S.Ns_nzocc(kpt),kpt)) + bsxfun(@times, -Chi_full * bsxfun(@times,integral_2',S.Atom(JJ_a).gamma_Jl) + ...
    	% 						                    DChi * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl) - DVeff_psi, transpose(S.occ(1:S.Ns_nzocc(kpt),kpt)));
    	% end
    
    	% for kpt = 1:S.tnkpt
    	% 	S = occupation_pert(S,scf_natmc);
    	% 	Beta_ph_diag = S.occ(1:S.Ns_nzocc(kpt),kpt) + 0.5*S.bet * S.occ(1:S.Ns_nzocc(kpt),kpt) .* (1-S.occ(1:S.Ns_nzocc(kpt),kpt)) .* (S.DlambdaF_dR(scf_natmc)./S.Dlambda_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc) - 1);
    
    	% 	rhs_strnhimr(:,1:S.Ns_nzocc(kpt),kpt) = rhs_strnhimr(:,1:S.Ns_nzocc(kpt),kpt) + bsxfun(@times,S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc), transpose(S.Dlambda_dR(1:S.Ns_nzocc(kpt),kpt,scf_natmc).*Beta_ph_diag));
    
    	% 	kpt_vec = S.kptgrid(kpt,:);
	        
	        %------------------------------Method-3 ends---------------------------------------
	        
	        
	        [DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kpt_vec);
	        
	        % Solve Sternheimer
	        % change the sign of Dpsi whenever psi changes sign for the guess
	        for Ns = 1:S.Ns_nzocc(kpt)
	            Hfun =  @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),kpt,kpt_vec,scf_natmc,fac,S,x,S.psi(:,1:S.Ns_nzocc(kpt),kpt));
	            %S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc) = gmres(Hfun,rhs_strnhimr(:,Ns,kpt),[],1e-6,3000);
	            S.Dpsi_dR_kq(:,Ns,kpt) = bicgstab(Hfun,rhs_strnhimr(:,Ns),S.linSolv_tol,3000,S.LapPreconL,S.LapPreconU,S.Dpsi_dR_kq(:,Ns,kpt));
	            %S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc) = aar_phonon(S,scf_natmc,Ns,kpt,kpt_vec,DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),rhs_strnhimr(:,Ns,kpt),S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc),S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc),S.linSolv_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,1);
	        end
	    end
	else

	    for kpt = 1:S.tnkpt
	        %-----------------------------------------Method-2---------------------------------------------------------
	        
	        %--------------------------------------------kq-----------------------------------------------
	        psi_kq = determine_psikq(S,kpt,1);
	        PsiW_k = transpose(bsxfun(@times, S.psi(:,1:S.Ns_nzocc(kpt),kpt), S.W));
	        PsiW_kq = transpose(bsxfun(@times, conj(psi_kq(:,1:S.Ns_nzocc_max+2)), S.W));
	        DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi(:,1:S.Ns_nzocc(kpt),kpt));
	        psi_DVeff_psi = PsiW_kq * DVeff_psi;
	        
	        % non-local psp (can be done only once during a scf)
	        fac = 1.0i;
	        JJ_a = ceil(scf_natmc/3);
	        Chi_full_k = zeros(S.N,S.Atom(JJ_a).angnum);
	        Chi_full_kq = zeros(S.N,S.Atom(JJ_a).angnum);
	        for img = 1:S.Atom(JJ_a).n_image_rc
	            img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	            phase_fac_k = exp(-1*dot(S.kptgrid(kpt,:),img_disp*fac));
	            phase_fac_kq = exp(-1*dot(S.kqptgrid(kpt,:),img_disp*fac));
	            Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_k;
	            Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_kq;
	        end
	        
	        

	        if Dir == 3
	        	DChi_k = blochGradient(S,S.kptgrid(kpt,:),Dir) * Chi_full_k;
	        	DChi_kq = blochGradient(S,S.kqptgrid(kpt,:),Dir) * Chi_full_kq;
	        elseif Dir == 1
	        	DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
	        	DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
	        	DChi_x_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),1) * Chi_full_kq;
	        	DChi_y_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),2) * Chi_full_kq;
	         	DChi_k = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_kq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			else
				DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
	        	DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
	        	DChi_x_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),1) * Chi_full_kq;
	        	DChi_y_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),2) * Chi_full_kq;
	         	DChi_k = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_kq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			end
	        
	        psikq_star_dchi = -PsiW_kq * DChi_kq;
	        chi_star_psik = PsiW_k * conj(Chi_full_k);
	        
	        psikq_star_chi = PsiW_kq * Chi_full_kq;
	        dchi_star_psik = -PsiW_k * conj(DChi_k);
	        
	        %temp =  integral_2 * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
	        psi_DVnl_psi = psikq_star_dchi * bsxfun(@times,transpose(chi_star_psik),S.Atom(JJ_a).gamma_Jl) + ...
	                       psikq_star_chi * bsxfun(@times,transpose(dchi_star_psik),S.Atom(JJ_a).gamma_Jl);

	        %S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt);
	        
	        psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        rhs_strnhimr(:,1:S.Ns_nzocc(kpt)) = psi_kq(:,1:S.Ns_nzocc_max+2) * (psi_DH_psi .* S.coeff_Pv_psi_kq(1:S.Ns_nzocc_max+2,1:S.Ns_nzocc(kpt),kpt)) - Chi_full_kq * bsxfun(@times,transpose(dchi_star_psik),S.Atom(JJ_a).gamma_Jl) + ...
	                                            DChi_kq * bsxfun(@times,transpose(chi_star_psik),S.Atom(JJ_a).gamma_Jl) - DVeff_psi;
	        
	        
	        [DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,S.kqptgrid(kpt,:));
	        
	        % Solve Sternheimer
	        % change the sign of Dpsi whenever psi changes sign for the guess
	        for Ns = 1:S.Ns_nzocc(kpt)
	            Hfun =  @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),kpt,S.kqptgrid(kpt,:),scf_natmc,fac,S,x,psi_kq(:,1:S.Ns_nzocc_max+2));
	            % S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc) = gmres(Hfun,rhs_strnhimr(:,Ns,kpt),[],1e-4,1000);
	            S.Dpsi_dR_kq(:,Ns,kpt) = bicgstab(Hfun,rhs_strnhimr(:,Ns),S.linSolv_tol,3000,S.LapPreconL,S.LapPreconU,S.Dpsi_dR_kq(:,Ns,kpt));
	            %S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc) = aar_phonon(S,scf_natmc,Ns,kpt,S.kqptgrid(kpt,:),DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),rhs_strnhimr(:,Ns,kpt),S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc),S.psi_kq(:,:,kpt),S.linSolv_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,1);
	        end
	    end
	    
	    %---------------------------------------minus kq-----------------------------------------------
	    for kpt = 1:S.tnkpt
	    	psi_mkq = determine_psikq(S,kpt,2);
	        PsiW_mk = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
	        PsiW_mkq = transpose(bsxfun(@times, conj(psi_mkq(:,1:S.Ns_nzocc_max+2)), S.W));
	        DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)));
	        psi_DVeff_psi = PsiW_mkq * DVeff_psi;
	        
	        % non-local psp (can be done only once during a scf)
	        fac = 1.0i;
	        JJ_a = ceil(scf_natmc/3);
	        Chi_full_mk = zeros(S.N,S.Atom(JJ_a).angnum);
	        Chi_full_mkq = zeros(S.N,S.Atom(JJ_a).angnum);
	        for img = 1:S.Atom(JJ_a).n_image_rc
	            img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	            phase_fac_mk = exp(-1*dot(-S.kptgrid(kpt,:),img_disp*fac));
	            phase_fac_mkq = exp(-1*dot(S.mkqptgrid(kpt,:),img_disp*fac));
	            Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mk;
	            Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mkq;
	        end

	        if Dir == 3
	        	DChi_mk = blochGradient(S,-S.kptgrid(kpt,:),Dir) * Chi_full_mk;
	        	DChi_mkq = blochGradient(S,S.mkqptgrid(kpt,:),Dir) * Chi_full_mkq;
	        elseif Dir == 1
	        	DChi_x_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),1) * Chi_full_mk;
	        	DChi_y_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),2) * Chi_full_mk;
	        	DChi_x_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),1) * Chi_full_mkq;
	        	DChi_y_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),2) * Chi_full_mkq;
	         	DChi_mk = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_mkq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			else
				DChi_x_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),1) * Chi_full_mk;
	        	DChi_y_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),2) * Chi_full_mk;
	        	DChi_x_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),1) * Chi_full_mkq;
	        	DChi_y_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),2) * Chi_full_mkq;
	         	DChi_mk = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_mkq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			end
	        
	        psimkq_star_dchi = -PsiW_mkq * DChi_mkq;
	        chi_star_psimk = PsiW_mk * conj(Chi_full_mk);
	        
	        psimkq_star_chi = PsiW_mkq * Chi_full_mkq;
	        dchi_star_psimk = -PsiW_mk * conj(DChi_mk);
	        
	        psi_DVnl_psi = psimkq_star_dchi * bsxfun(@times,transpose(chi_star_psimk),S.Atom(JJ_a).gamma_Jl) + ...
	                       psimkq_star_chi * bsxfun(@times,transpose(dchi_star_psimk),S.Atom(JJ_a).gamma_Jl);
	        
	        
	        psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        rhs_strnhimr(:,1:S.Ns_nzocc(kpt)) = psi_mkq(:,1:S.Ns_nzocc_max+2) * (psi_DH_psi .* S.coeff_Pv_psi_mkq(1:S.Ns_nzocc_max+2,1:S.Ns_nzocc(kpt),kpt)) - Chi_full_mkq * bsxfun(@times,transpose(dchi_star_psimk),S.Atom(JJ_a).gamma_Jl) + ...
	                                            DChi_mkq * bsxfun(@times,transpose(chi_star_psimk),S.Atom(JJ_a).gamma_Jl) - DVeff_psi;
	        
	        [DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,S.mkqptgrid(kpt,:));
	        
	        % Solve Sternheimer
	        % change the sign of Dpsi whenever psi changes sign for the guess
	        for Ns = 1:S.Ns_nzocc(kpt)
	            Hfun =  @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),kpt,S.mkqptgrid(kpt,:),scf_natmc,fac,S,x,psi_mkq(:,1:S.Ns_nzocc_max+2));
	            % S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc) = gmres(Hfun,rhs_strnhimr(:,Ns,kpt),[],1e-4,1000);
	            S.Dpsi_dR_mkq(:,Ns,kpt) = bicgstab(Hfun,rhs_strnhimr(:,Ns),S.linSolv_tol,3000,S.LapPreconL,S.LapPreconU,S.Dpsi_dR_mkq(:,Ns,kpt));
	            %S.Dpsi_dR_mkq(:,Ns,kpt,scf_natmc) = aar_phonon(S,scf_natmc,Ns,kpt,S.mkqptgrid(kpt,:),DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),rhs_strnhimr(:,Ns,kpt),conj(S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc)),S.psi_mkq(:,:,kpt),S.linSolv_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,1);
	        end
	    end
	    
	end

end


function S = sternheimer_kparal(S,scf_natmc,qpt)

	Dpsi_dR = zeros(S.N,S.Ns_nzocc_max,S.tnkpt);
	Dlambda_dR = zeros(S.Ns_nzocc_max,S.tnkpt);

	Dir = rem(scf_natmc,3);
	if Dir == 0
	    Dir = 3;
	end

	if S.Kdeltaq0 == 1
		
		LASTN = maxNumCompThreads(1);
		parfor (kpt = 1:S.tnkpt, S.num_worker_heuristic)
			rhs_strnhimr = zeros(S.N,S.Ns_nzocc(kpt));
			
	        
	        %-----------------------------------Method-2-------------------------------------------------
	        %============================================================================================
	        PsiW = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
	        DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi(:,1:S.Ns_nzocc(kpt),kpt));
	        psi_DVeff_psi = PsiW * DVeff_psi;
	        
	        % non-local psp (can be done only once during a scf)
	        kpt_vec = S.kptgrid(kpt,:);
	        fac = 1.0i;
	        JJ_a = ceil(scf_natmc/3);
	        Chi_full = zeros(S.N,S.Atom(JJ_a).angnum);
	        for img = 1:S.Atom(JJ_a).n_image_rc
	            img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	            phase_fac = exp(-1*dot(kpt_vec,img_disp*fac));
	            Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac;
	        end
	        
	        if Dir == 3
	        	DChi = blochGradient(S,kpt_vec,Dir) * Chi_full;
	        elseif Dir == 1
	        	DChi_x_temp = blochGradient(S,kpt_vec,1)*Chi_full;
	 			DChi_y_temp = blochGradient(S,kpt_vec,2)*Chi_full;
	         	DChi = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			else
				DChi_x_temp = blochGradient(S,kpt_vec,1)*Chi_full;
	 			DChi_y_temp = blochGradient(S,kpt_vec,2)*Chi_full;
	         	DChi = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			end
	        
	        integral_1 = conj(PsiW) * conj(Chi_full);
	        integral_2 = -1*PsiW*DChi;
	        
	        temp =  integral_2 * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
	        psi_DVnl_psi = temp + temp';
	       
	        %S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt);
	        
	        psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        %Dlambda_dR(:,kpt) = diag(psi_DH_psi);
	        coeff_Pv_Psi = S.coeff_Pv_psi_kq(1:S.Ns_nzocc(kpt),1:S.Ns_nzocc(kpt),kpt);% + diag(ones(S.Ns_nzocc(kpt),1));
	        rhs_strnhimr(:,:) = S.psi(:,1:S.Ns_nzocc(kpt),kpt) * (psi_DH_psi .* coeff_Pv_Psi) - Chi_full * bsxfun(@times,integral_2',S.Atom(JJ_a).gamma_Jl) + ...
	                                                DChi * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl) - DVeff_psi;


	        %---------------------------Method-2 ends-----------------------------------------------------        
	        [DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kpt_vec);
	        
	        % Solve Sternheimer
	        % change the sign of Dpsi whenever psi changes sign for the guess
	        [Dpsi_dR(:,:,kpt),Dlambda_dR(:,kpt)] = linSolve(S,DL11,DL22,DL33,DG1,DG2,DG3,kpt,kpt_vec,scf_natmc,fac,S.psi(:,1:S.Ns_nzocc(kpt),kpt),rhs_strnhimr(:,:),Dpsi_dR(:,:,kpt),S.Dpsi_dR_kq(:,:,kpt),Dlambda_dR(:,kpt),psi_DH_psi);
	        % for Ns = 1:S.Ns_nzocc_max
	        % 	Hfun =  @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),kpt,kpt_vec,scf_natmc,fac,S,x,S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc));
	        %     % Dpsi_dR(:,Ns,kpt) = gmres(Hfun,rhs_strnhimr(:,Ns,kpt),[],1e-4,1000);
	        %     Dpsi_dR(:,:,kpt) =  pcg(Hfun,rhs_strnhimr(:,Ns,kpt),S.linSolv_tol,3000,S.LapPreconL,S.LapPreconU,S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc));
	        %     %Dpsi_dR(:,Ns,kpt) = aar_phonon(S,scf_natmc,Ns,kpt,kpt_vec,DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),rhs_strnhimr(:,Ns,kpt),S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc),S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc),S.linSolv_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,1);
	        % end
	    end
	    maxNumCompThreads(LASTN);
	    S.Dlambda_dR(:,:,scf_natmc) = Dlambda_dR;
	    S.Dpsi_dR_kq(:,:,:) = Dpsi_dR;

	else

	    LASTN = maxNumCompThreads(1);
		parfor (kpt = 1:S.tnkpt, S.num_worker_heuristic)
			rhs_strnhimr = zeros(S.N,S.Ns_nzocc(kpt));
	        %-----------------------------------------Method-2---------------------------------------------------------
	        
	        %--------------------------------------------kq-----------------------------------------------
	        psi_kq = determine_psikq(S,kpt,1);
	        PsiW_k = transpose(bsxfun(@times, S.psi(:,1:S.Ns_nzocc(kpt),kpt), S.W));
	        PsiW_kq = transpose(bsxfun(@times, conj(psi_kq(:,1:S.Ns_nzocc_max+2)), S.W));
	        DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),S.psi(:,1:S.Ns_nzocc(kpt),kpt));
	        psi_DVeff_psi = PsiW_kq * DVeff_psi;
	        
	        % non-local psp (can be done only once during a scf)
	        fac = 1.0i;
	        JJ_a = ceil(scf_natmc/3);
	        Chi_full_k = zeros(S.N,S.Atom(JJ_a).angnum);
	        Chi_full_kq = zeros(S.N,S.Atom(JJ_a).angnum);
	        for img = 1:S.Atom(JJ_a).n_image_rc
	            img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	            phase_fac_k = exp(-1*dot(S.kptgrid(kpt,:),img_disp*fac));
	            phase_fac_kq = exp(-1*dot(S.kqptgrid(kpt,:),img_disp*fac));
	            Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_k;
	            Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_kq;
	        end
	        
	       if Dir == 3
	        	DChi_k = blochGradient(S,S.kptgrid(kpt,:),Dir) * Chi_full_k;
	        	DChi_kq = blochGradient(S,S.kqptgrid(kpt,:),Dir) * Chi_full_kq;
	        elseif Dir == 1
	        	DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
	        	DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
	        	DChi_x_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),1) * Chi_full_kq;
	        	DChi_y_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),2) * Chi_full_kq;
	         	DChi_k = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_kq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			else
				DChi_x_temp_k = blochGradient(S,S.kptgrid(kpt,:),1) * Chi_full_k;
	        	DChi_y_temp_k = blochGradient(S,S.kptgrid(kpt,:),2) * Chi_full_k;
	        	DChi_x_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),1) * Chi_full_kq;
	        	DChi_y_temp_kq = blochGradient(S,S.kqptgrid(kpt,:),2) * Chi_full_kq;
	         	DChi_k = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_kq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_k(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_kq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			end
	        
	        psikq_star_dchi = -PsiW_kq * DChi_kq;
	        chi_star_psik = PsiW_k * conj(Chi_full_k);
	        
	        psikq_star_chi = PsiW_kq * Chi_full_kq;
	        dchi_star_psik = -PsiW_k * conj(DChi_k);
	        
	        %temp =  integral_2 * bsxfun(@times,transpose(integral_1),S.Atom(JJ_a).gamma_Jl);
	        psi_DVnl_psi = psikq_star_dchi * bsxfun(@times,transpose(chi_star_psik),S.Atom(JJ_a).gamma_Jl) + ...
	                       psikq_star_chi * bsxfun(@times,transpose(dchi_star_psik),S.Atom(JJ_a).gamma_Jl);

	        %S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) = S.psi(:,1:S.Ns_nzocc(kpt),kpt);
	        
	        psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        rhs_strnhimr(:,:) = psi_kq(:,1:S.Ns_nzocc_max+2) * (psi_DH_psi .* S.coeff_Pv_psi_kq(1:S.Ns_nzocc_max+2,1:S.Ns_nzocc(kpt),kpt)) - Chi_full_kq * bsxfun(@times,transpose(dchi_star_psik),S.Atom(JJ_a).gamma_Jl) + ...
	                            DChi_kq * bsxfun(@times,transpose(chi_star_psik),S.Atom(JJ_a).gamma_Jl) - DVeff_psi;
	        
	        
	        [DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,S.kqptgrid(kpt,:));
	        
	        % Solve Sternheimer
	        % change the sign of Dpsi whenever psi changes sign for the guess
	        [Dpsi_dR(:,:,kpt),~] = linSolve(S,DL11,DL22,DL33,DG1,DG2,DG3,kpt,S.kqptgrid(kpt,:),scf_natmc,fac,psi_kq(:,1:S.Ns_nzocc_max+2),rhs_strnhimr(:,:),Dpsi_dR(:,:,kpt),S.Dpsi_dR_kq(:,:,kpt),Dlambda_dR(:,kpt),psi_DH_psi);
	        
	        % for Ns = 1:S.Ns_nzocc(kpt)
	        %     Hfun =  @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),kpt,S.kqptgrid(kpt,:),scf_natmc,fac,S,x,S.psi_kq(:,:,kpt));
	        %     % Dpsi_dR(:,Ns,kpt) = gmres(Hfun,rhs_strnhimr(:,Ns,kpt),[],1e-4,1000);
	        %     Dpsi_dR(:,Ns,kpt) = pcg(Hfun,rhs_strnhimr(:,Ns,kpt),S.linSolv_tol,3000,S.LapPreconL,S.LapPreconU,S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc));
	        %     %Dpsi_dR(:,Ns,kpt) = aar_phonon(S,scf_natmc,Ns,kpt,S.kqptgrid(kpt,:),DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),rhs_strnhimr(:,Ns,kpt),S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc),S.psi_kq(:,:,kpt),S.linSolv_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,1);
	        % end
	    end
	    maxNumCompThreads(LASTN);
		S.Dpsi_dR_kq(:,:,:) = Dpsi_dR;
	    
	    %---------------------------------------minus kq-----------------------------------------------
	    LASTN = maxNumCompThreads(1);
		parfor (kpt = 1:S.tnkpt, S.num_worker_heuristic)
			rhs_strnhimr = zeros(S.N,S.Ns_nzocc(kpt));
			psi_mkq = determine_psikq(S,kpt,2);

	        PsiW_mk = transpose(bsxfun(@times, conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)), S.W));
	        PsiW_mkq = transpose(bsxfun(@times, conj(psi_mkq(:,1:S.Ns_nzocc_max+2)), S.W));
	        DVeff_psi = bsxfun(@times,S.Dveff_dR(:,scf_natmc),conj(S.psi(:,1:S.Ns_nzocc(kpt),kpt)));
	        psi_DVeff_psi = PsiW_mkq * DVeff_psi;
	        
	        % non-local psp (can be done only once during a scf)
	        fac = 1.0i;
	        JJ_a = ceil(scf_natmc/3);
	        Chi_full_mk = zeros(S.N,S.Atom(JJ_a).angnum);
	        Chi_full_mkq = zeros(S.N,S.Atom(JJ_a).angnum);
	        for img = 1:S.Atom(JJ_a).n_image_rc
	            img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	            phase_fac_mk = exp(-1*dot(-S.kptgrid(kpt,:),img_disp*fac));
	            phase_fac_mkq = exp(-1*dot(S.mkqptgrid(kpt,:),img_disp*fac));
	            Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mk;
	            Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Chi_full_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * phase_fac_mkq;
	        end
	        
	        if Dir == 3
	        	DChi_mk = blochGradient(S,-S.kptgrid(kpt,:),Dir) * Chi_full_mk;
	        	DChi_mkq = blochGradient(S,S.mkqptgrid(kpt,:),Dir) * Chi_full_mkq;
	        elseif Dir == 1
	        	DChi_x_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),1) * Chi_full_mk;
	        	DChi_y_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),2) * Chi_full_mk;
	        	DChi_x_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),1) * Chi_full_mkq;
	        	DChi_y_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),2) * Chi_full_mkq;
	         	DChi_mk = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_mkq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,1)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,1)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			else
				DChi_x_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),1) * Chi_full_mk;
	        	DChi_y_temp_mk = blochGradient(S,-S.kptgrid(kpt,:),2) * Chi_full_mk;
	        	DChi_x_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),1) * Chi_full_mkq;
	        	DChi_y_temp_mkq = blochGradient(S,S.mkqptgrid(kpt,:),2) * Chi_full_mkq;
	         	DChi_mk = zeros(S.N,S.Atom(JJ_a).angnum);
	         	DChi_mkq = zeros(S.N,S.Atom(JJ_a).angnum);

				for img = 1:S.Atom(JJ_a).n_image_rc
					img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
					fac1 = -img_disp(2)/S.L2;
					fac2 = -img_disp(3)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mk(S.Atom(JJ_a).rcImage(img).rc_pos,:);
					DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) = DChi_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(1,2)*DChi_x_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:) + ROT(2,2)*DChi_y_temp_mkq(S.Atom(JJ_a).rcImage(img).rc_pos,:);
				end
			end
	        
	        psimkq_star_dchi = -PsiW_mkq * DChi_mkq;
	        chi_star_psimk = PsiW_mk * conj(Chi_full_mk);
	        
	        psimkq_star_chi = PsiW_mkq * Chi_full_mkq;
	        dchi_star_psimk = -PsiW_mk * conj(DChi_mk);
	        
	        psi_DVnl_psi = psimkq_star_dchi * bsxfun(@times,transpose(chi_star_psimk),S.Atom(JJ_a).gamma_Jl) + ...
	                       psimkq_star_chi * bsxfun(@times,transpose(dchi_star_psimk),S.Atom(JJ_a).gamma_Jl);
	        
	        
	        psi_DH_psi = psi_DVeff_psi + psi_DVnl_psi;
	        rhs_strnhimr(:,:) = psi_mkq(:,1:S.Ns_nzocc_max+2) * (psi_DH_psi .* S.coeff_Pv_psi_mkq(1:S.Ns_nzocc_max+2,1:S.Ns_nzocc(kpt),kpt)) - Chi_full_mkq * bsxfun(@times,transpose(dchi_star_psimk),S.Atom(JJ_a).gamma_Jl) + ...
	                            DChi_mkq * bsxfun(@times,transpose(chi_star_psimk),S.Atom(JJ_a).gamma_Jl) - DVeff_psi;
	        
	        [DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,S.mkqptgrid(kpt,:));
	        
	        % Solve Sternheimer
	        % change the sign of Dpsi whenever psi changes sign for the guess
	        [Dpsi_dR(:,:,kpt),~] = linSolve(S,DL11,DL22,DL33,DG1,DG2,DG3,kpt,S.mkqptgrid(kpt,:),scf_natmc,fac,psi_mkq(:,1:S.Ns_nzocc_max+2),rhs_strnhimr(:,:),Dpsi_dR(:,:,kpt),S.Dpsi_dR_mkq(:,:,kpt),Dlambda_dR(:,kpt),psi_DH_psi);
	        % for Ns = 1:S.Ns_nzocc(kpt)
	        %     Hfun =  @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),kpt,S.mkqptgrid(kpt,:),scf_natmc,fac,S,x,S.psi_mkq(:,:,kpt));
	        %     % Dpsi_dR_mkq(:,Ns,kpt) = gmres(Hfun,rhs_strnhimr(:,Ns,kpt),[],1e-4,1000);
	        %     Dpsi_dR_mkq(:,Ns,kpt) = pcg(Hfun,rhs_strnhimr(:,Ns,kpt),S.linSolv_tol,3000,S.LapPreconL,S.LapPreconU,S.Dpsi_dR_mkq(:,Ns,kpt,scf_natmc));
	        %     %Dpsi_dR_mkq(:,Ns,kpt) = aar_phonon(S,scf_natmc,Ns,kpt,S.mkqptgrid(kpt,:),DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),rhs_strnhimr(:,Ns,kpt),conj(S.Dpsi_dR_kq(:,Ns,kpt,scf_natmc)),S.psi_mkq(:,:,kpt),S.linSolv_tol,1000,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU,1);
	        % end
	    end
	    maxNumCompThreads(LASTN);
	    S.Dpsi_dR_mkq(:,:,:) = Dpsi_dR;
	    
	end

end


function [Dpsi_dR,Dlambda_dR] = linSolve(S,DL11,DL22,DL33,DG1,DG2,DG3,kpt,kpt_vec,scf_natmc,fac,psi_v,rhs_strnhimr,Dpsi_dR,Dpsi_dR_guess,Dlambda_dR,psi_DH_psi);
	for Ns = 1:S.Ns_nzocc(kpt)
		Hfun =  @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,S.EigVal(Ns,kpt),kpt,kpt_vec,scf_natmc,fac,S,x,psi_v);
		Dpsi_dR(:,Ns) =  bicgstab(Hfun,rhs_strnhimr(:,Ns),S.linSolv_tol,3000,S.LapPreconL,S.LapPreconU,Dpsi_dR_guess(:,Ns));%gmres(Hfun,rhs_strnhimr(:,Ns),[],1e-6,3000);%
	end
	Dlambda_dR(1:S.Ns_nzocc(kpt)) = diag(psi_DH_psi);
end



function [x] = aar_phonon(S,scf_natmc,Ns,kpt,kptvec,DL11,DL22,DL33,DG1,DG2,DG3,eigval,b,x_prev,psi_v,tol,max_iter,omega,beta,m,p,L,U,opt)

	% 	fac = 1.0i;
	

	% Hnlx = -0.5*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,x_prev,S)) + bsxfun(@times,S.Veff-eigval,x_prev) + S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc) * ( (S.psi_n(:,1:S.Ns_nzocc(kpt),kpt,scf_natmc)' * (x_prev .* S.W)) ); %.* (abs(S.EigVal(1:S.Ns_nzocc(kpt),kpt) - eigval) < 1e-3)
	% for JJ_a = 1:S.n_atm
	% 	Chi_X_mult = zeros(S.Atom(JJ_a).angnum,1);
	% 	for img = 1:S.Atom(JJ_a).n_image_rc
	% 		img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	% 		phase_fac = exp(dot(kptvec,img_disp*fac));
	% 		Chi_X_mult = Chi_X_mult + ( bsxfun(@times, S.Atom(JJ_a).rcImage(img).Chi_mat, S.W(S.Atom(JJ_a).rcImage(img).rc_pos)) )' * x_prev(S.Atom(JJ_a).rcImage(img).rc_pos) * phase_fac;
	% 	end
	% 	% Chi_X_mult = Chi_X_mult .* S.Atom(J).gamma_Jl;
	% 	Chi_X_mult = bsxfun(@times,Chi_X_mult, S.Atom(JJ_a).gamma_Jl);
	% 	for img = 1:S.Atom(JJ_a).n_image_rc
	% 		img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
	% 		phase_fac = exp(-1*dot(kptvec,img_disp*fac));
	% 		Hnlx(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Hnlx(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * Chi_X_mult * phase_fac;
	% 	end
	% end



	% AAR begins here


	N=length(b);
	DX=zeros(N,m);
	DF=zeros(N,m);
	nb = 1/norm(b);
	% if isempty(x_guess)
	% 	x_guess = ones(N,1);  % initial guess vector
	% end
	% x_prev = x_guess;
	relres = 1+tol;
	count = 1;


	if opt == 1
		fac = 1.0i;
		Ax = @(x) H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,eigval,kpt,kptvec,scf_natmc,fac,S,x,psi_v);
	else
		Ax = @(x) Lap_times_x(DL11,DL22,DL33,DG1,DG2,DG3,x,S);
	end

	%Hnlx = Ax(x_prev);

	while count<=max_iter && relres > tol

		Hnlx = Ax(x_prev);

		res1 = (b-Hnlx) ;
		relres = norm(res1)*nb;
		res = U\(L\(res1)); 
		
		% STORE HISTORY
		if count>1
			DX(:,mod(count-2,m)+1) = x_prev-Xold;
			DF(:,mod(count-2,m)+1) = res-Fold;
		end
		Xold = x_prev;
		Fold = res;
		
		% UPDATE ITERATE
		if rem(count,p)~=0   % RICHARDSON UPDATE  
			x_new = x_prev + omega*res; 
		else % ANDERSON UPDATE, apply every p iters
			x_new = x_prev + beta*res - (DX + beta*DF)*(pinv(DF'*DF)*(DF'*res));
		end
		
		x = x_prev;
		x_prev = x_new;
		count = count + 1;
	end

	if count-1 == max_iter
		fprintf('Sternheimer: AAR exceeded maximum iterations and converged to a relative residual of %g. \n',relres);
	else
		fprintf('Sternheimer: AAR converged to a relative residual of %g in %d iterations.\n',relres,count-1);
	end

end
	

function [Hnlx] = H_times_x(DL11,DL22,DL33,DG1,DG2,DG3,eigval,kpt,kptvec,scf_natmc,fac,S,x_prev,psi_v)
	Hnlx = -0.5*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,x_prev,S)) + bsxfun(@times,S.Veff-eigval,x_prev) + psi_v * (psi_v' * bsxfun(@times,x_prev,S.W));
	
	for JJ_a = 1:S.n_atm
		Chi_X_mult = zeros(S.Atom(JJ_a).angnum,1);
		for img = 1:S.Atom(JJ_a).n_image_rc
			img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
			phase_fac = exp(dot(kptvec,img_disp*fac));
			Chi_X_mult = Chi_X_mult + ( bsxfun(@times, S.Atom(JJ_a).rcImage(img).Chi_mat, S.W(S.Atom(JJ_a).rcImage(img).rc_pos)) )' * x_prev(S.Atom(JJ_a).rcImage(img).rc_pos) * phase_fac;
		end
		% Chi_X_mult = Chi_X_mult .* S.Atom(J).gamma_Jl;
		Chi_X_mult = bsxfun(@times,Chi_X_mult, S.Atom(JJ_a).gamma_Jl);
		for img = 1:S.Atom(JJ_a).n_image_rc
			img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
			phase_fac = exp(-1*dot(kptvec,img_disp*fac));
			Hnlx(S.Atom(JJ_a).rcImage(img).rc_pos,:) = Hnlx(S.Atom(JJ_a).rcImage(img).rc_pos,:) + S.Atom(JJ_a).rcImage(img).Chi_mat * Chi_X_mult * phase_fac;
		end
	end
end

function [Lapx] = Lap_times_x(DL11,DL22,DL33,DG1,DG2,DG3,x,S)
	Lapx = lapVec(DL11,DL22,DL33,DG1,DG2,DG3,x,S);
end


function S = electrostatics_local(S,qpt)

	qpt_vec = S.qptgrid(qpt,:);
	fac = 1.0i;
	
	M_dyn_temp1 = ones(3*S.n_atm,S.N);
	M_dyn_temp2 = ones(S.N,3*S.n_atm);
	Vec_temp = zeros(S.N,1);
	SC_row_x = zeros(1,S.N);
	SC_row_y = zeros(1,S.N);
	SC_row_z = zeros(1,S.N);
	SC_col_x = zeros(S.N,1);
	SC_col_y = zeros(S.N,1);
	SC_col_z = zeros(S.N,1);
	M_mass = zeros(3*S.n_atm,1);


	Dphi_x = S.grad_1*(S.phi);
	Dphi_y = S.grad_2*(S.phi);
	Dphi_z = S.grad_3*(S.phi);

	DVc_x = S.grad_1*(S.V_c);
	DVc_y = S.grad_2*(S.V_c);
	DVc_z = S.grad_3*(S.V_c);

	Db_plus_bref_x = S.grad_1*(S.b + S.b_ref);
	Db_plus_bref_y = S.grad_2*(S.b + S.b_ref);
	Db_plus_bref_z = S.grad_3*(S.b + S.b_ref);

	% Store derivatives needed for electrostatics
	%=============================================
	count_typ = 1;
	count_typ_atms = 1;
	for JJ_a = 1:S.n_atm % loop over all the atoms

		% Create mass matrix
		M_mass(3*(JJ_a-1)+1:3*(JJ_a-1)+3,1) = S.Atm(count_typ).Mass;

		SC_row_x = 0 * SC_row_x;
		SC_row_y = 0 * SC_row_y;
		SC_row_z = 0 * SC_row_z;
		SC_col_x = 0 * SC_col_x;
		SC_col_y = 0 * SC_col_y;
		SC_col_z = 0 * SC_col_z;
		
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
		xx_img = (-n_image_xl : n_image_xr) * S.L1 + x0;
		yy_img = (-n_image_yl : n_image_yr) * S.L2 + y0;
		zz_img = (-n_image_zl : n_image_zr) * S.L3 + z0;
		[XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = ndgrid(xx_img,yy_img,zz_img);

		% Loop over all image(s) of atom JJ_a (including atom JJ_a)
		for count_image = 1:n_image_total

			Vec_temp = 0 * Vec_temp;
		
			% Atom position of the image
			x0_i = XX_IMG_3D(count_image);
			y0_i = YY_IMG_3D(count_image);
			z0_i = ZZ_IMG_3D(count_image);

			img_disp = [x0_i-x0 y0_i-y0 z0_i-z0];
 			phase_fac = exp(dot(qpt_vec,img_disp*fac));
		
			% Starting and ending indices of b-region
			ii_s = ceil (((x0_i-S.xin) - S.Atm(count_typ).rb_x)/S.dx) + 1;
			ii_e = floor(((x0_i-S.xin) + S.Atm(count_typ).rb_x)/S.dx) + 1;
			jj_s = ceil (((y0_i-S.yin) - S.Atm(count_typ).rb_y)/S.dy) + 1;
			jj_e = floor(((y0_i-S.yin) + S.Atm(count_typ).rb_y)/S.dy) + 1;
			kk_s = ceil (((z0_i-S.zin) - S.Atm(count_typ).rb_z)/S.dz) + 1;
			kk_e = floor(((z0_i-S.zin) + S.Atm(count_typ).rb_z)/S.dz) + 1;
		
		
			 % Check if the b-region is inside the domain in Dirichlet BC
			 % direction
			isInside = (S.BCx == 0 || (S.BCx == 1 && (ii_s>=1) && (ii_e<=S.Nx))) && ...
				       (S.BCy == 0 || (S.BCy == 1 && (jj_s>=1) && (jj_e<=S.Ny))) && ...
				       (S.BCz == 0 || (S.BCz == 1 && (kk_s>=1) && (kk_e<=S.Nz)));
			% assert(isInside,'Error: Atom too close to boundary for b calculation');
			if ~isInside
				fprintf(' WARNING: Atom %d too close to boundary for b calculation\n',JJ_a);
			end
			ii_s = max(ii_s,1);
			ii_e = min(ii_e,S.Nx);
			jj_s = max(jj_s,1);
			jj_e = min(jj_e,S.Ny);
			kk_s = max(kk_s,1);
			kk_e = min(kk_e,S.Nz);

			% xx = S.xin + (ii_s-3*S.FDn-1:ii_e+3*S.FDn-1)*S.dx;
			% yy = S.yin + (jj_s-3*S.FDn-1:jj_e+3*S.FDn-1)*S.dy;
			% zz = S.zin + (kk_s-3*S.FDn-1:kk_e+3*S.FDn-1)*S.dz;

			xx = S.xin + (ii_s-2*S.FDn-1:ii_e+2*S.FDn-1)*S.dx;
			yy = S.yin + (jj_s-2*S.FDn-1:jj_e+2*S.FDn-1)*S.dy;
			zz = S.zin + (kk_s-2*S.FDn-1:kk_e+2*S.FDn-1)*S.dz;
			
			[XX_3D,YY_3D,ZZ_3D] = ndgrid(xx,yy,zz);
		
			% Find distances
			dd = calculateDistance(XX_3D,YY_3D,ZZ_3D,x0_i,y0_i,z0_i,S);
		
			% Pseudopotential at grid points through interpolation
			V_PS = zeros(size(dd));
			IsLargeThanRmax = dd > S.Atm(count_typ).r_grid_vloc(end);
			V_PS(IsLargeThanRmax) = -S.Atm(count_typ).Z;
			V_PS(~IsLargeThanRmax) = interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).r_grid_vloc.*S.Atm(count_typ).Vloc, dd(~IsLargeThanRmax), 'spline');
			
			V_PS = V_PS./dd;
			V_PS(dd<S.Atm(count_typ).r_grid_vloc(2)) = S.Atm(count_typ).Vloc(1);
			
			% Reference potential at grid points
			rc_ref = S.rc_ref; % WARNING: Might need smaller if pseudocharges overlap
			V_PS_ref = zeros(size(dd));
			I_ref = dd<rc_ref;
			V_PS_ref(~I_ref) = -(S.Atm(count_typ).Z)./dd(~I_ref);
			V_PS_ref(I_ref) = -S.Atm(count_typ).Z*(9*dd(I_ref).^7-30*rc_ref*dd(I_ref).^6 ...
				+28*rc_ref*rc_ref*dd(I_ref).^5-14*(rc_ref^5)*dd(I_ref).^2+12*rc_ref^7)/(5*rc_ref^8);

			% Isolated atom electron density at grid points through interpolation
			rho_isolated_atom = interp1(S.Atm(count_typ).r_grid_rho, S.Atm(count_typ).rho_isolated_guess, dd, 'spline');
			rho_isolated_atom(dd > S.Atm(count_typ).r_grid_rho(end)) = 0; 
		
			% Pseudocharge density
			II = 1+S.FDn : size(V_PS,1)-S.FDn;
			JJ = 1+S.FDn : size(V_PS,2)-S.FDn;
			KK = 1+S.FDn : size(V_PS,3)-S.FDn;
			
			% Calculate bJ and bJ_ref
			bJ = pseudochargeDensity_atom(V_PS,II,JJ,KK,xx(1),S);
			bJ_ref = pseudochargeDensity_atom(V_PS_ref,II,JJ,KK,xx(1),S);
			
			bJ = (-1/(4*pi))*bJ;
			bJ_ref = (-1/(4*pi))*bJ_ref;

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% Calculate the gradient of pseudocharges
			II = 1+2*S.FDn : size(V_PS,1)-2*S.FDn;
			JJ = 1+2*S.FDn : size(V_PS,2)-2*S.FDn;
			KK = 1+2*S.FDn : size(V_PS,3)-2*S.FDn;
			[dbJ_x, dbJ_y, dbJ_z] = dpseudopot(bJ,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
			[dbJ_ref_x, dbJ_ref_y, dbJ_ref_z] = dpseudopot(bJ_ref,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
			[dVJ_x, dVJ_y, dVJ_z] = dpseudopot(V_PS,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
			[dVJ_ref_x, dVJ_ref_y, dVJ_ref_z] = dpseudopot(V_PS_ref,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
			[drho_x, drho_y, drho_z] = dpseudopot(rho_isolated_atom,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);

			dbJ_x_qpt = dbJ_x * phase_fac; dbJ_y_qpt = dbJ_y * phase_fac; dbJ_z_qpt = dbJ_z * phase_fac;
			dbJ_ref_x_qpt = dbJ_ref_x * phase_fac; dbJ_ref_y_qpt = dbJ_ref_y * phase_fac; dbJ_ref_z_qpt = dbJ_ref_z * phase_fac;
			dVJ_x_qpt = dVJ_x * phase_fac; dVJ_y_qpt =  dVJ_y * phase_fac; dVJ_z_qpt = dVJ_z * phase_fac;
			dVJ_ref_x_qpt = dVJ_ref_x * phase_fac; dVJ_ref_y_qpt = dVJ_ref_y * phase_fac; dVJ_ref_z_qpt = dVJ_ref_z * phase_fac;
			drho_x_qpt = drho_x * phase_fac; drho_y_qpt = drho_y * phase_fac; drho_z_qpt = drho_z * phase_fac;

			[II_rb,JJ_rb,KK_rb] = ndgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
			Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;
			M_dyn_elec_temp = zeros(3,3);

			M_dyn_elec_temp(1,1) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_x(II,JJ,KK).*(Dphi_x(Rowcount_rb)))));
			M_dyn_elec_temp(1,2) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_y(II,JJ,KK).*(Dphi_x(Rowcount_rb)))));
			M_dyn_elec_temp(1,3) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_z(II,JJ,KK).*(Dphi_x(Rowcount_rb)))));
			M_dyn_elec_temp(2,1) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_x(II,JJ,KK).*(Dphi_y(Rowcount_rb)))));
			M_dyn_elec_temp(2,2) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_y(II,JJ,KK).*(Dphi_y(Rowcount_rb)))));
			M_dyn_elec_temp(2,3) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_z(II,JJ,KK).*(Dphi_y(Rowcount_rb)))));
			M_dyn_elec_temp(3,1) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_x(II,JJ,KK).*(Dphi_z(Rowcount_rb)))));
			M_dyn_elec_temp(3,2) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_y(II,JJ,KK).*(Dphi_z(Rowcount_rb)))));
			M_dyn_elec_temp(3,3) = - sum(sum(sum(S.W(Rowcount_rb).*dbJ_z(II,JJ,KK).*(Dphi_z(Rowcount_rb)))));

			M_dyn_elec_temp(1,1) = M_dyn_elec_temp(1,1) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_x(II,JJ,KK)+dbJ_ref_x(II,JJ,KK)).*(DVc_x(Rowcount_rb)))));
			M_dyn_elec_temp(1,2) = M_dyn_elec_temp(1,2) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_y(II,JJ,KK)+dbJ_ref_y(II,JJ,KK)).*(DVc_x(Rowcount_rb)))));
			M_dyn_elec_temp(1,3) = M_dyn_elec_temp(1,3) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_z(II,JJ,KK)+dbJ_ref_z(II,JJ,KK)).*(DVc_x(Rowcount_rb)))));
			M_dyn_elec_temp(2,1) = M_dyn_elec_temp(2,1) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_x(II,JJ,KK)+dbJ_ref_x(II,JJ,KK)).*(DVc_y(Rowcount_rb)))));
			M_dyn_elec_temp(2,2) = M_dyn_elec_temp(2,2) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_y(II,JJ,KK)+dbJ_ref_y(II,JJ,KK)).*(DVc_y(Rowcount_rb)))));
			M_dyn_elec_temp(2,3) = M_dyn_elec_temp(2,3) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_z(II,JJ,KK)+dbJ_ref_z(II,JJ,KK)).*(DVc_y(Rowcount_rb)))));
			M_dyn_elec_temp(3,1) = M_dyn_elec_temp(3,1) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_x(II,JJ,KK)+dbJ_ref_x(II,JJ,KK)).*(DVc_z(Rowcount_rb)))));
			M_dyn_elec_temp(3,2) = M_dyn_elec_temp(3,2) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_y(II,JJ,KK)+dbJ_ref_y(II,JJ,KK)).*(DVc_z(Rowcount_rb)))));
			M_dyn_elec_temp(3,3) = M_dyn_elec_temp(3,3) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(dbJ_z(II,JJ,KK)+dbJ_ref_z(II,JJ,KK)).*(DVc_z(Rowcount_rb)))));

			M_dyn_elec_temp(1,1) = M_dyn_elec_temp(1,1) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_x(II,JJ,KK)+dVJ_ref_x(II,JJ,KK)).*(Db_plus_bref_x(Rowcount_rb)))));
			M_dyn_elec_temp(1,2) = M_dyn_elec_temp(1,2) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_y(II,JJ,KK)+dVJ_ref_y(II,JJ,KK)).*(Db_plus_bref_x(Rowcount_rb)))));
			M_dyn_elec_temp(1,3) = M_dyn_elec_temp(1,3) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_z(II,JJ,KK)+dVJ_ref_z(II,JJ,KK)).*(Db_plus_bref_x(Rowcount_rb)))));
			M_dyn_elec_temp(2,1) = M_dyn_elec_temp(2,1) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_x(II,JJ,KK)+dVJ_ref_x(II,JJ,KK)).*(Db_plus_bref_y(Rowcount_rb)))));
			M_dyn_elec_temp(2,2) = M_dyn_elec_temp(2,2) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_y(II,JJ,KK)+dVJ_ref_y(II,JJ,KK)).*(Db_plus_bref_y(Rowcount_rb)))));
			M_dyn_elec_temp(2,3) = M_dyn_elec_temp(2,3) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_z(II,JJ,KK)+dVJ_ref_z(II,JJ,KK)).*(Db_plus_bref_y(Rowcount_rb)))));
			M_dyn_elec_temp(3,1) = M_dyn_elec_temp(3,1) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_x(II,JJ,KK)+dVJ_ref_x(II,JJ,KK)).*(Db_plus_bref_z(Rowcount_rb)))));
			M_dyn_elec_temp(3,2) = M_dyn_elec_temp(3,2) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_y(II,JJ,KK)+dVJ_ref_y(II,JJ,KK)).*(Db_plus_bref_z(Rowcount_rb)))));
			M_dyn_elec_temp(3,3) = M_dyn_elec_temp(3,3) - 0.5*sum(sum(sum(S.W(Rowcount_rb).*(-dVJ_z(II,JJ,KK)+dVJ_ref_z(II,JJ,KK)).*(Db_plus_bref_z(Rowcount_rb)))));

			
			fac1 = -(y0-y0_i)/S.L2;
			fac2 = -(z0-z0_i)/S.L3;
			ROT = (S.RotM1^fac1) * (S.RotM2^fac2); 
			S.M_dyn_elec(3*(JJ_a-1)+1:3*(JJ_a-1)+3,3*(JJ_a-1)+1:3*(JJ_a-1)+3) = S.M_dyn_elec(3*(JJ_a-1)+1:3*(JJ_a-1)+3,3*(JJ_a-1)+1:3*(JJ_a-1)+3) + ROT' * M_dyn_elec_temp * ROT;
			

			Vec_temp(Rowcount_rb) = dbJ_x_qpt(II,JJ,KK); Vec_temp_x = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = dbJ_y_qpt(II,JJ,KK); Vec_temp_y = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = dbJ_z_qpt(II,JJ,KK); Vec_temp_z = reshape(Vec_temp,[],1);
			P = [Vec_temp_x Vec_temp_y Vec_temp_z] * ROT;
			Vec_temp_x = P(:,1); Vec_temp_y = P(:,2); Vec_temp_z = P(:,3);

			S.Db_dR(:,3*(JJ_a-1)+1) = S.Db_dR(:,3*(JJ_a-1)+1) - Vec_temp_x;
			S.Db_dR(:,3*(JJ_a-1)+2) = S.Db_dR(:,3*(JJ_a-1)+2) - Vec_temp_y;
			S.Db_dR(:,3*(JJ_a-1)+3) = S.Db_dR(:,3*(JJ_a-1)+3) - Vec_temp_z;

			Vec_temp(Rowcount_rb) = drho_x_qpt(II,JJ,KK); Vec_temp_x = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = drho_y_qpt(II,JJ,KK); Vec_temp_y = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = drho_z_qpt(II,JJ,KK); Vec_temp_z = reshape(Vec_temp,[],1);
			P = [Vec_temp_x Vec_temp_y Vec_temp_z] * ROT;
			Vec_temp_x = P(:,1); Vec_temp_y = P(:,2); Vec_temp_z = P(:,3);

			S.Drho_dR(:,3*(JJ_a-1)+1) = S.Drho_dR(:,3*(JJ_a-1)+1) - Vec_temp_x;
			S.Drho_dR(:,3*(JJ_a-1)+2) = S.Drho_dR(:,3*(JJ_a-1)+2) - Vec_temp_y;
			S.Drho_dR(:,3*(JJ_a-1)+3) = S.Drho_dR(:,3*(JJ_a-1)+3) - Vec_temp_z;

			Vec_temp(Rowcount_rb) = dbJ_x_qpt(II,JJ,KK) + dbJ_ref_x_qpt(II,JJ,KK); Vec_temp_x = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = dbJ_y_qpt(II,JJ,KK) + dbJ_ref_y_qpt(II,JJ,KK); Vec_temp_y = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = dbJ_z_qpt(II,JJ,KK) + dbJ_ref_z_qpt(II,JJ,KK); Vec_temp_z = reshape(Vec_temp,[],1);
			P = [Vec_temp_x Vec_temp_y Vec_temp_z] * ROT;
			Vec_temp_x = P(:,1); Vec_temp_y = P(:,2); Vec_temp_z = P(:,3);

			SC_row_x = SC_row_x + transpose(Vec_temp_x);
			SC_row_y = SC_row_y + transpose(Vec_temp_y);
			SC_row_z = SC_row_z + transpose(Vec_temp_z);

			Vec_temp(Rowcount_rb) = dVJ_ref_x_qpt(II,JJ,KK) - dVJ_x_qpt(II,JJ,KK); Vec_temp_x = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = dVJ_ref_y_qpt(II,JJ,KK) - dVJ_y_qpt(II,JJ,KK); Vec_temp_y = reshape(Vec_temp,[],1);
			Vec_temp(Rowcount_rb) = dVJ_ref_z_qpt(II,JJ,KK) - dVJ_z_qpt(II,JJ,KK); Vec_temp_z = reshape(Vec_temp,[],1);
			P = [Vec_temp_x Vec_temp_y Vec_temp_z] * ROT;
			Vec_temp_x = P(:,1); Vec_temp_y = P(:,2); Vec_temp_z = P(:,3);

			SC_col_x = SC_col_x + Vec_temp_x;
			SC_col_y = SC_col_y + Vec_temp_y;
			SC_col_z = SC_col_z + Vec_temp_z;
		end
		
		% Check if same type of atoms are over
		if count_typ_atms == S.Atm(count_typ).n_atm_typ
			count_typ_atms = 1;
			count_typ = count_typ + 1;
		else
			count_typ_atms = count_typ_atms + 1;
		end

		M_dyn_temp1(3*(JJ_a-1)+1,:) = 0.5*SC_row_x;
		M_dyn_temp1(3*(JJ_a-1)+2,:) = 0.5*SC_row_y;
		M_dyn_temp1(3*(JJ_a-1)+3,:) = 0.5*SC_row_z;

		M_dyn_temp2(:,3*(JJ_a-1)+1) = SC_col_x.*S.W;
		M_dyn_temp2(:,3*(JJ_a-1)+2) = SC_col_y.*S.W;
		M_dyn_temp2(:,3*(JJ_a-1)+3) = SC_col_z.*S.W;


	end % end of loop over atoms

	%l = S.M_dyn_elec

	temp = conj(M_dyn_temp1) * M_dyn_temp2;
	S.M_dyn_elec = S.M_dyn_elec + temp + temp';

	%l = S.M_dyn_elec
	%eig(S.M_dyn_elec)
	%aa
    
	S.M_mass = diag(M_mass);

end




function [psi_kq] = determine_psikq(S,kpt,typ)
    psi_kq = zeros(S.N,S.Nev);
    TOL = 1e-8;
    N_c = round(2*pi/S.L2);

    if typ == 1

		chkmatch1 = abs(S.kqptgrid(kpt,1) - S.kptgrid(:,1)) < TOL & abs(S.kqptgrid(kpt,2) - S.kptgrid(:,2)) < TOL & abs(S.kqptgrid(kpt,3) - S.kptgrid(:,3)) < TOL;
        chkmatch2 = abs(S.kqptgrid(kpt,1) + S.kptgrid(:,1)) < TOL & (abs(S.kqptgrid(kpt,2) + S.kptgrid(:,2) - N_c) < TOL | abs(S.kqptgrid(kpt,2) + S.kptgrid(:,2)) < TOL) & abs(S.kqptgrid(kpt,3) + S.kptgrid(:,3)) < TOL;
                
		if sum(chkmatch1) == 1
		    psi_kq(:,:) = S.psi(:,:,chkmatch1 == 1);
		elseif sum(chkmatch2) == 1
		    psi_kq(:,:) = conj(S.psi(:,:,chkmatch2 == 1));
		else
			error('wrong kpt');
		end

	elseif typ == 2

        chkmatch1 = abs(S.mkqptgrid(kpt,1) - S.kptgrid(:,1)) < TOL & abs(S.mkqptgrid(kpt,2) - S.kptgrid(:,2)) < TOL & abs(S.mkqptgrid(kpt,3) - S.kptgrid(:,3)) < TOL;
        chkmatch2 = abs(S.mkqptgrid(kpt,1) + S.kptgrid(:,1)) < TOL & (abs(S.mkqptgrid(kpt,2) + S.kptgrid(:,2) - N_c) < TOL | abs(S.mkqptgrid(kpt,2) + S.kptgrid(:,2)) < TOL) & abs(S.mkqptgrid(kpt,3) + S.kptgrid(:,3)) < TOL;
                
        if sum(chkmatch1) == 1
            psi_kq(:,:) = S.psi(:,:,chkmatch1 == 1);
        elseif sum(chkmatch2) == 1
            psi_kq(:,:) = conj(S.psi(:,:,chkmatch2 == 1));
        else
           error('wrong kpt');
        end
    end
end


function S = const_for_FFT(S)
	w2 = S.w2;
	FDn = S.FDn;
	N1 = S.Nx;
	N2 = S.Ny;
	N3 = S.Nz;
	dx = S.dx;
	dy = S.dy;
	dz = S.dz;

	[I,J,K] = meshgrid(0:(N1-1),0:(N2-1),0:(N3-1));

	dx2 = dx*dx; dy2 = dy*dy; dz2 = dz*dz;

	% alpha follows conjugate even space
	alpha = w2(1)*(1/dx2+1/dy2+1/dz2).*ones(N1,N2,N3);
	for k=1:FDn
	    alpha = alpha + w2(k+1)*2.*(cos(2*pi*I*k/N1)./dx2 + cos(2*pi*J*k/N2)./dy2 + cos(2*pi*K*k/N3)./dz2);
	end

	V = S.L1*S.L2*S.L3;
	R_c = (3*V/(4*pi))^(1/3);
	alpha(1,1,1) = -2/R_c^2;

	const = 1 - cos(R_c*sqrt(-1*alpha));
	const(1,1,1) = 1;
	S.const_by_alpha = const./alpha;
end










function [V] = poissonSolve_FFT(S,rhs)
	% if(S.BC ~= 2)
	%     error("Must use Periodic BC\n");
	% end

	% t1 = tic;
	f = -4 * pi * rhs;
	f = reshape(f,S.Nx,S.Ny,S.Nz);
	g_hat = fftn(f);
	V = ifftn(g_hat.*S.const_by_alpha);
	V = V(:);
	% fprintf(' Poisson problem solved by FFT took %fs\n',toc(t1));
end



function f = poisson_RHS(S,rhs)
	% @brief	Poisson_RHS evaluates the right hand side of the poisson equation, 
	%           including the boundary condtions, i.e. f = -4 * pi * ( rho + b - d) 
	%           for cluster system, while for periodic system it's just 
	%           f = -4 * pi * (rho + b).
	f = -4 * pi * rhs;

	% for charged systems, add a uniform background charge so that total charge
	% is 0
	if S.NetCharge ~= 0
		unif_bkgd_chrg = S.NetCharge / (S.L1*S.L2*S.L3);
		f = f + (-4 * pi) * unif_bkgd_chrg;
		% fprintf(2,'total charge + background charge = %d\n',dot(S.W,rho+S.b+unif_bkgd_chrg));
	end

	if(S.BC == 1)
		% For cluster systems, we need to include boundary conditions d
		% RR = S.RR_AUG(S.isIn);
		for l = 0:S.l_cut
			multipole_moment(l+1).Qlm = sum(repmat(S.RR.^l .* (rhs) .* S.W,1,2*l+1).* S.SH(l+1).Ylm )';
		end

		% Calculate phi using multipole expansion
		phi = zeros(size(S.RR_AUG_3D));
		for l = 0 : S.l_cut
			denom = (2*l+1)*S.RR_AUG_3D.^(l+1);
			for m = -l : l
				Ylm_AUG_3D = reshape(S.SH(l+1).Ylm_AUG(:,m+l+1),size(phi));
				phi = phi + Ylm_AUG_3D .* multipole_moment(l+1).Qlm(m+l+1) ./ denom;
			end
		end
		phi = 4 * pi * phi;
		phi(S.isIn) = 0;
	elseif S.BC == 2
		return;
	elseif S.BC == 3
		% find boundary conditions for periodic in 2D, dirichlet in 1D
		% Calculate phi
	    phi = zeros(S.Nx+2*S.FDn, S.Ny+2*S.FDn, S.Nz+2*S.FDn);

		cellsize  = [S.L1, S.L2, S.L3];
		gridsizes = [S.Nx, S.Ny, S.Nz];
		meshsizes = [S.dx, S.dy, S.dz];
		bcs       = [S.BCx,S.BCy,S.BCz];

		% find which direction has Dirichlet BC
		dir_Z = find(bcs == 1); 
		dir_X = mod(dir_Z, 3) + 1;
		dir_Y = mod(dir_X, 3) + 1;

		% once we find the direction, we assume that direction is the Z
		% direction, the other two directions are then called X, Y
		NX = gridsizes(dir_X); NY = gridsizes(dir_Y); NZ = gridsizes(dir_Z);
		LX = cellsize(dir_X);  LY = cellsize(dir_Y);  LZ = cellsize(dir_Z);
		dX = meshsizes(dir_X); dY = meshsizes(dir_Y); dZ = meshsizes(dir_Z);
		A_XY = LX * LY; % area of (x',y') surface

		% reshape rho to 3D 
		rho = reshape(rhs, S.Nx, S.Ny, S.Nz); % note here after rho = rho + b

		% permute rho so that the new z' direction has Dirichlet BC
		new_order = [dir_X, dir_Y, dir_Z]; % a permutation of [1,2,3]
		[~, reverse_order] = sort(new_order); % reverse order to get back
		rho = permute(rho,new_order);
		phi = permute(phi,new_order);
		% rho = permute(rho,reverse_order); % this will recover original

		% find fft of rho in the new X, Y directions
		% rho is a 3D matrix, this applies 2D fft to each dim higher than 2
		rho_hat = fft2(rho); 
		rho_hat(1,1,:) = 0; % remove the 0 frequency term (constant term in real-space)

		% rho_zp_bar = zeros(Nzp,1);
		% sum over x' and y' directions, \int (\rho) dxdy
		%rho_zp_bar = sum(sum(rho, dir_xp), dir_yp) * dxp * dyp; 
		rho_Z_av = sum(sum(rho, 1), 2) * dX * dY; 
		rho_Z_av = rho_Z_av(:);
		Z = (0 : NZ-1) * dZ;
		Z_bc = [-S.FDn:-1,NZ:NZ+S.FDn-1] * dZ;

		% a matrix of size Nz' x Nz', a(i,j) =  z_i - z_j
		%ZZp = abs(Z_bc.' - Z); 
		ZZp = abs(colminusrow(Z_bc.', Z)); 

		% V_av(Z) = int (rho_Z_av * |Z - Z'|) dZ', note this can be simplified
		% for charge neutrual systems to int (rho_Z_av * Z') dZ', indep of Z
		V_Z_av = (-2*pi/A_XY*dZ) * sum(bsxfun(@times, rho_Z_av', ZZp),2);

		% wave vectors
		GX = ifftshift( -floor(NX/2) : ceil(NX/2)-1 ) * (2*pi/NX);
		GY = ifftshift( -floor(NY/2) : ceil(NY/2)-1 ) * (2*pi/NY);
		[GX3D, GY3D, Z3D] = ndgrid(GX, GY, Z);
		GR3D = sqrt(GX3D.^2 + GY3D.^2);
		V_hat = zeros(NX,NY,2*S.FDn);
		for k = 1:S.ksymlength(Z_bc)  % in fact, no need to calculate for all Z values
			V_hat(:,:,k) = V_hat(:,:,k) + (2*pi*dZ) * ...
				sum(rho_hat ./ GR3D .* exp(-abs(Z3D - Z_bc(k)) .* GR3D) , 3);
		end
		V_hat(1,1,:) = 0; % remove the 0 frequency term (constant term in real-space)
		
		II = (1+S.FDn):(NX+S.FDn);
		JJ = (1+S.FDn):(NY+S.FDn);

		% put phi back to the 3D matrix
		ind_z = [1:S.FDn,NZ+S.FDn+1:NZ+2*S.FDn];
		phi(II,JJ,ind_z) = reshape(kron(V_Z_av', ones(NX,NY)),NX,NY,S.ksymlength(Z_bc));

		% remove 0 before ifft2 to include more terms
		phi(II,JJ,ind_z) = phi(II,JJ,ind_z) + 0*ifft2(V_hat); 
		% phi(S.isIn) = 0; % any value inside the fd-grid will be set to 0

		% permute phi back to original directions
		phi = real(phi);
		phi = permute(phi,reverse_order);
	elseif S.BC == 4
		% find boundary conditions for periodic in 1D, dirichlet in 2D
		% Calculate phi
		phi = zeros(S.Nx+2*S.FDn, S.Ny+2*S.FDn, S.Nz+2*S.FDn);

		cellsize  = [S.L1, S.L2, S.L3];
		gridsizes = [S.Nx, S.Ny, S.Nz];
		meshsizes = [S.dx, S.dy, S.dz];
		bcs       = [S.BCx,S.BCy,S.BCz];

		% find which direction has Periodic BC
		dir_Z = find(bcs == 0); 
		dir_X = mod(dir_Z, 3) + 1;
		dir_Y = mod(dir_X, 3) + 1;

		% once we find the direction, we assume that direction is the Z
		% direction, the other two directions are then called X, Y
		NX = gridsizes(dir_X); NY = gridsizes(dir_Y); NZ = gridsizes(dir_Z);
		LX = cellsize(dir_X);  LY = cellsize(dir_Y);  LZ = cellsize(dir_Z);
		dX = meshsizes(dir_X); dY = meshsizes(dir_Y); dZ = meshsizes(dir_Z);
		% A_XY = LX * LY; % area of (x',y') surface

		% reshape rho to 3D 
		rho = reshape(rhs, S.Nx, S.Ny, S.Nz); % note here after rho = rho + b

		% permute rho so that the new Z direction has Periodic BC
		new_order = [dir_X, dir_Y, dir_Z]; % a permutation of [1,2,3]
		[~, reverse_order] = sort(new_order); % reverse order to get back
		rho = permute(rho,new_order);
		phi = permute(phi,new_order);
		% rho = permute(rho,reverse_order); % this will recover original

		% sum over Z direction, \int (\rho) dz / LZ
		rho_XY_av = sum(rho,3) * (dZ/LZ); % NX * NY
		rho_XY_av = rho_XY_av(:);
		X = (0 : NX-1) * dX;
		Y = (0 : NY-1) * dY;

		[XX,YY] = ndgrid(X,Y);

		I_ex = 1:NX+2*S.FDn;
		J_ex = 1:NY+2*S.FDn;
		[II_ex,JJ_ex] = ndgrid(I_ex,J_ex);

		% flag for points inside the domain
		% isIn = ones(NX+2*S.FDn,NY+2*S.FDn);
		% isIn(1:S.FDn,:) = 0;
		% isIn(:,1:S.FDn) = 0;
		% isIn(S.FDn+NX+1:end,:) = 0;
		% isIn(:,S.FDn+NY+1:end) = 0;
		isIn = zeros(NX+2*S.FDn,NY+2*S.FDn);
		isIn(S.FDn+1:NX+S.FDn, S.FDn+1:NY+S.FDn) = 1;

		% positions of nodes outside the domain
		isbc = find(~isIn);

		XX_bc = (II_ex(isbc)-1-S.FDn) * dX;
		YY_bc = (JJ_ex(isbc)-1-S.FDn) * dY;

		% ln((X-X')^2 + (Y-Y')^2)
		% ln_RRp = log((XX_bc - XX(:)').^2 + (YY_bc - YY(:)').^2 );
		ln_RRp = log(colminusrow(XX_bc,XX(:)').^2 + colminusrow(YY_bc,YY(:)').^2);
		
		% V_XY_av = int (-ln((X-X')^2 + (Y-Y')^2) * rho_XY_av) dX'dY'
		V_XY_av = (-dX*dY) * sum(bsxfun(@times, rho_XY_av', ln_RRp),2);

		V_XY_av_full = zeros(NX+2*S.FDn,NY+2*S.FDn);
		V_XY_av_full(isbc) = V_XY_av;

		phi = bsxfun(@plus, phi, V_XY_av_full);

		% permute phi back to original directions
		phi = real(phi);
		phi = permute(phi,reverse_order);
	end
end

% tool function colminusrow
function xmy = colminusrow(x,y)
% A column vector x minus a row vector.
% In Matlab versions after R2018b, it's just x - y
if (size(x,2) ~= 1)
	error('Error: the first vector must be a column vector');
end
if (size(y,1) ~= 1)
	error('Error: the second vector must be a row vector');
end

[xx,yy] = ndgrid(x,y);
xmy = xx - yy;

end



function [eigVal, eigVec] = eigSolve_k(S,kptvec)

	kpt0 = 1;
	mind = 10;
	for kpt = 1:S.tnkpt
		dif = abs(S.kptgrid(kpt,1) - kptvec(1)) + abs(S.kptgrid(kpt,2) - kptvec(2)) + abs(S.kptgrid(kpt,3) - kptvec(3));
		if (dif < mind)
			mind = dif;
			kpt0 = kpt;
		end
	end

	% Heff = spdiags(S.Veff(:,spin),0,S.N,S.N);
	Heff = S.Veff;
    rng('default'); % Initialize random number generator
    eigVec = S.psi(:,:,kpt0);
    upper_bound_guess_vecs = zeros(S.N,1);
    eigVal = zeros(S.Nev,1);
	%opts = struct('maxit', 10000, 'tol', 1e-6, 'p', S.Nev+10, 'v0', rand(S.N,1), 'isreal', true);
	opts = struct('maxit', 1000, 'tol', 1e-5, 'v0', rand(S.N,1));
	[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kptvec);
	Hfun = @(x) h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,x,S,kptvec);
	if ~(isreal(DL11) && isreal(DL22) && isreal(DL33))
		opts.isreal = false;
	end

	% Upper bound estimator
	% [eigVec(:,1:S.Ns_nzocc_max),eigsm,flag] = eigs(Hfun,S.N,S.Ns_nzocc_max,'sr',opts);
	% eigVal(1:S.Ns_nzocc_max,1) = diag(eigsm);
	% flag

	% scfac = 1 ./ sqrt(sum(repmat(S.W,1,S.Nev) .* (eigVec .* conj(eigVec)),1));
	% eigVec = bsxfun(@times, eigVec, scfac);
	
	[upper_bound_guess_vecs, bup] = (eigs(Hfun,S.N,1,'lr',opts));
	bup = real(bup) + opts.tol;
	% Lower bound estimator
	a0 = real(eigs(Hfun,S.N,1,'sr',opts)) - opts.tol;
	lambda_cutoff = max(S.EigVal(:,kpt0)) + 0.15;
	
	% fprintf('filter cutoff = %f, lower bound = %f, upper bound = %f\n',lambda_cutoff(ks),a0(ks),bup(ks));
	%count = 0;
	err = 1;
	eigVal_prev = S.EigVal(:,kpt0);
	while (err > 1e-8)
		% Chebyshev filtering
		eigVec = chebyshev_filter(eigVec,S.npl,lambda_cutoff,bup,a0,DL11,DL22,DL33,DG1,DG2,DG3,Heff,S,kptvec);
		eigVec = orth(eigVec);
		Nev1 = size(eigVec,2);    % WARNING: ORTH(psi) might change the size of psi, that's why we update Nev
		assert(Nev1 == S.Nev,'Number of states have changed');

		% Subspace Hamiltonian
		Hs = eigVec' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Heff,eigVec,S,kptvec);

		% Solve subspace eigenproblem,
		if S.cell_typ < 3
			Hs = 0.5 * (Hs + Hs');
		end
		[Q, Q1] = eig(Hs);
		eigVal = real(diag(Q1)); % WARNING: Taking real part only!

		% subspace rotation
		eigVec = eigVec * Q;

		% Normalize psi, s.t. integral(psi_new' * psi_new) = 1
		scfac = 1 ./ sqrt(sum(repmat(S.W,1,S.Nev) .* (eigVec .* conj(eigVec)),1));
		% psi(:,:,ks) = psi(:,:,ks) * diag(scfac);
		eigVec = bsxfun(@times, eigVec, scfac);

		err = max(abs(eigVal(1:S.Ns_nzocc_max,1) - eigVal_prev(1:S.Ns_nzocc_max,1)));


		eigVal_prev = eigVal;
	end

end







