function [atompos_next,AtmForce,S] = electronicGroundState(atom_pos,S)
% @brief    electronicGroundState(atom_pos,S) calculates the electronic ground
%           state atomic force corresponding to the given atom positions.
%
% @param atom_pos       Given atom positions (stretched into a (3*n_atoms)-by-1 vector)
% @param atompos_next   Atom position for the next relaxation/MD step (in cartesian coordinates)
% @param AtmForce       Atomic force as a 3*N x 1 vector 
% @authors	Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%==============================================================================================

tic_relax = tic;
fprintf('\n');
fprintf(' ###############################################################\n');
fprintf(' Relaxation step number: %d \n', S.Relax_iter);

if (S.RelaxFlag == 1 || S.MDFlag)
	% Reshaping atom positions
	atom_pos = transpose(reshape(atom_pos,3,[]));
	% update un-mapped atomic positions, for charge extrapolation
	if S.cell_typ <= 2
	    S.atom_pos_tp1 = S.atom_pos_tp1 + (atom_pos - S.Atoms_old);
	else
	    S.atom_pos_tp1 = atom_pos;
	end    
	% S.dV_tp1 = S.dV;
	% Convert atom coordinates from cartesian to cell coordinates
	S.Atoms = coordinateTransformation(S, atom_pos, 'cart2noncart_dis');
end



% perform charge extrapolation (before rho_at is updated)
if (S.RelaxFlag == 1 || S.MDFlag)
    % Check position of atom near the boundary and apply wraparound in case of PBC
    S = check_atomlocation(S);
	S = chargeExtrapolation(S);
end

% Pseudocharge (and reference), sum atomic charge density, self energy 
% (and reference), electrostatic correction 
% S.b,S.b_ref,S.Eself,S.Eself_ref,S.rho_at,S.E_corr,S.V_c, S.NegCharge,
% S.PosCharge, S.NetCharge
t_calc_b = tic;

S = calculate_b_guessRho_Eself(S);

fprintf(' Time for b calculation: %.3f seconds.\n',toc(t_calc_b));

% set up guess electron density (guess rho)
S = initElectrondensity(S);

% Calculate nonlocal projectors	
S.Atom = calculate_nloc_projector(S);

% Self-consistent Field (SCF) method
S = scf(S);
	
%save('rho.mat','-struct','S','rho');
S.S_Debug.relax(S.Relax_iter).Eself = S.Eself;
S.S_Debug.relax(S.Relax_iter).Eself_ref = S.Eself_ref;
S.S_Debug.relax(S.Relax_iter).E_corr = S.E_corr;

if abs(1-S.occ(1))>1e-6 || abs(S.occ(end))>1e-6
	fprintf('[\b Warning: No. of states is not enough!]\b \n');
	S.S_Debug.relax(S.Relax_iter).occ_check = 1; % 1 means not satisfied
else
	S.S_Debug.relax(S.Relax_iter).occ_check = 0;
end
	
% Etotal = evaluateTotalEnergy(EigVal,occ,rho,S.b,phi,Vxc,S.W,S.bet,S.Eself,S.E_corr,1) ;

fprintf('\n');
fprintf(' **********************************************************\n');
fprintf(' *          Energy per unit cell = %11.9f Ha.       *\n', S.Etotal);
fprintf(' *          Energy per atom = %11.9f Ha.            *\n', S.Etotal / S.n_atm);
fprintf(' **********************************************************\n');

% write to output file
outfname = S.outfname;
fileID = fopen(outfname,'a');
if (fileID == -1) 
	error('\n Cannot open file "%s"\n',outfname);
end 
fprintf(fileID,'====================================================================\n');
fprintf(fileID,'                                Energy                              \n');
fprintf(fileID,'====================================================================\n');
fprintf(fileID,'Free energy per atom               :%18.10E (Ha/atom)\n', S.Etotal / S.n_atm);
fprintf(fileID,'Total free energy                  :%18.10E (Ha)\n', S.Etotal);
fprintf(fileID,'Band structure energy              :%18.10E (Ha)\n', S.Eband);
fprintf(fileID,'Exchange correlation energy        :%18.10E (Ha)\n', S.Exc);
fprintf(fileID,'Self and correction energy         :%18.10E (Ha)\n', S.E_corr-S.Eself);
fprintf(fileID,'Entropy*kb*T                       :%18.10E (Ha)\n', S.Eent);
fprintf(fileID,'Fermi level                        :%18.10E (Ha)\n', S.lambda_f);
if S.nspin ~= 1
	fprintf(fileID,'Net Magnetization                  :%18.10E \n', S.netM);
end
fclose(fileID);


% Atomic force calculation
tic_force = tic;
S.force = atomicForce(S);
%S.force = zeros(S.n_atm,3);
force_mat = S.force;
sz_fmat = size(force_mat);
force_corr = sum(force_mat, 1) / sz_fmat(1);
if S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5
	force_corr(1) = 0.0; % x - component
	force_corr(2) = 0.0; % y - component
end
S.force = force_mat - repmat(force_corr, sz_fmat(1),1);

% Apply constraint on atoms by making forces on fixed atoms zero
if S.RelaxFlag || S.MDFlag
	S.force = S.force .* S.mvAtmConstraint;
end

fprintf(' ***********************************************************\n');
fprintf(' *                      Atomic Force                       *\n');
fprintf(' ***********************************************************\n');
fprintf(' Drift free forces (Ha/Bohr):\n');
S.S_Debug.relax(S.Relax_iter).force = S.force;
disp(S.force);
S.abs_force = sqrt(sum(abs(S.force).^2,2));
fprintf(' Max magnitude of forces (Ha/Bohr):');
disp(max(S.abs_force));

t_force = toc(tic_force);
fprintf('\n Time for calculating forces: %f s.\n', t_force);

% write forces into .static file if required
if (S.PrintForceFlag == 1 && S.MDFlag == 0 && S.RelaxFlag == 0) 
	staticfname = S.staticfname;
	fileID = fopen(staticfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',staticfname);
	end
	fprintf(fileID, 'Atomic forces (Ha/Bohr):\n');
	fprintf(fileID, '%18.10f %18.10f %18.10f\n', S.force');
	fclose(fileID);
end

% write force magnitude and timing info. to .out file
outfname = S.outfname;
fileID = fopen(outfname,'a');
if (fileID == -1) 
	error('\n Cannot open file "%s"\n',outfname);
end

force_magnitude = sqrt(sum(S.force .* S.force, 2));
avgF = sum(force_magnitude) / length(force_magnitude);
maxF = max(force_magnitude);

fprintf(fileID,'Average force                      :%18.10E (Ha/Bohr)\n',avgF);
fprintf(fileID,'Maximum force                      :%18.10E (Ha/Bohr)\n',maxF);
fprintf(fileID,'Time for force calculation         :  %.3f (sec)\n',t_force);
fclose(fileID);


% Assign output of function
%AtmForce = reshape(S.force',[],1);
%AtmForce = - AtmForce;


% Perform stress and pressure calculation
if S.Calc_stress
	t1 = tic;
	S.Stress = evaluateStress(S);
    if (S.BC == 2)
        S.Stress = S.Stress * S.Ha_Bohr3_GPa;
        S.Pressure = -trace(S.Stress)/3.0;
        fprintf('\n[\b"Stress in GPa"\n\n\n]\b');
        disp(S.Stress);
        S.stress_dimred = eye(3);
    elseif (S.BC == 3)
        if(S.BCx == 0 && S.BCy == 0)
            %S.Stress = [sigma(1,1) sigma(1,2); sigma(2,1) sigma(2,2)];
            L_vac = S.L3;
            S.stress_dimred = [1 0 0; 0 1 0];
        elseif(S.BCx == 0 && S.BCz == 0)
            %S.Stress = [sigma(1,1) sigma(1,3); sigma(3,1) sigma(3,3)];
            L_vac = S.L2;
            S.stress_dimred = [1 0 0; 0 0 1];
        else
            %S.Stress = [sigma(2,2) sigma(2,3); sigma(3,2) sigma(3,3)];
            L_vac = S.L1;
            S.stress_dimred = [0 1 0; 0 0 1];
        end
        
        S.Stress = S.stress_dimred * S.Stress * S.stress_dimred';
        fprintf('\n[\b"Stress in Ha/Bohr**2"\n\n\n]\b');
        disp(S.Stress);
        
        stress_inGPa = (S.Stress/L_vac) * S.Ha_Bohr3_GPa;
        fprintf('\n[\b"Stress in GPa (equivalent to all periodic system)"\n\n\n]\b');
        disp(stress_inGPa);
    elseif (S.BC == 4)
        if(S.BCx == 0)
            %S.Stress = sigma(1,1);
            A_vac = S.L2*S.L3;
            S.stress_dimred = [1 0 0];
        elseif(S.BCy == 0)
            %S.Stress = sigma(2,2);
            A_vac = S.L1*S.L3;
            S.stress_dimred = [0 1 0];
        else
            %S.Stress = sigma(3,3);
            A_vac = S.L1*S.L2;
            S.stress_dimred = [0 0 1];
        end
        
        S.Stress = S.stress_dimred * S.Stress * S.stress_dimred';
        fprintf('\n[\b"Stress in Ha/Bohr"\n\n\n]\b');
        disp(S.Stress);
        
        stress_inGPa = (S.Stress/A_vac) * S.Ha_Bohr3_GPa;
        fprintf('\n[\b"Stress in GPa (equivalent to all periodic system)"\n\n\n]\b');
        disp(stress_inGPa);
    elseif (S.BC >= 5 && S.BC <= 7)
    	fprintf('\n[\b"Stress in Ha/Bohr"\n\n\n]\b');
        disp(S.Stress);
        A_vac = pi * (S.xout^2 - S.xin^2);
        stress_inGPa = (S.Stress/A_vac) * S.Ha_Bohr3_GPa;
        fprintf('\n[\b"Stress in GPa (equivalent to all periodic system)"\n\n\n]\b');
        disp(stress_inGPa);
    end
    
    fprintf('\n Time for calculating stress: %f s.\n\n', toc(t1));

	outfname = S.outfname;
	fileID = fopen(outfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',outfname);
	end
    if (S.MDFlag == 0 && S.RelaxFlag == 0)
        if (S.BC == 2)
            fprintf(fileID, 'Stress (GPa)                       :\n');
            fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.Stress(1,:));
            fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.Stress(2,:));
            fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.Stress(3,:));
        elseif (S.BC == 3)
            fprintf(fileID, 'Stress (Ha/Bohr**2)                :\n');
            fprintf(fileID,'%18.10E %18.10E \n',S.Stress(1,:));
            fprintf(fileID,'%18.10E %18.10E \n',S.Stress(2,:));
            
            fprintf(fileID, 'Stress equiv. all periodic (GPa)   :\n');
            fprintf(fileID,'%18.10E %18.10E \n',stress_inGPa(1,:));
            fprintf(fileID,'%18.10E %18.10E \n',stress_inGPa(2,:));
        elseif (S.BC >= 4 && S.BC <= 7)
            fprintf(fileID, 'Stress (Ha/Bohr)                   :\n');
            fprintf(fileID,'%18.10E \n',S.Stress(1,:));
            
            fprintf(fileID, 'Stress equiv. all periodic (GPa)   :\n');
            fprintf(fileID,'%18.10E \n',stress_inGPa(1,:));
        end
    end
    if (S.BC == 2)
        fprintf(fileID,'Pressure                           :%18.10E (GPa)\n',S.Pressure);
        fprintf(fileID,'Max stress                         :%18.10E (GPa)\n',max(max(abs(S.Stress))));
    elseif (S.BC == 3)
        fprintf(fileID,'Max stress                         :%18.10E (Ha/Bohr**2)\n',max(max(abs(S.Stress))));
        fprintf(fileID,'Max stress equiv. all periodic     :%18.10E (GPa)\n',max(max(abs(stress_inGPa))));
    elseif (S.BC >= 4 && S.BC <= 7)
        fprintf(fileID,'Max stress                         :%18.10E (Ha/Bohr)\n',max(abs(S.Stress)));
        fprintf(fileID,'Max stress equiv. all periodic     :%18.10E (GPa)\n',max(abs(stress_inGPa)));
    end
	fprintf(fileID,'Time for stress calculation        :  %.3f (sec)\n',toc(t1));
	
	fclose(fileID);        
end


if (S.Calc_pres && ~S.Calc_stress)
	t1 = tic;
	Pressure = evaluatePressure(S);
	PinGPa = Pressure*S.Ha_Bohr3_GPa;
	S.Pressure = PinGPa;  % store pressure
	fprintf('\n[\b"Pressure %f Ha/Bohr^3 %f GPa"\n\n\n]\b',Pressure,PinGPa);
	fprintf('\n Time for calculating pressure: %f s.\n\n', toc(t1));

	outfname = S.outfname;
	fileID = fopen(outfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',outfname);
	end
	fprintf(fileID,'Pressure                           :%18.10E (GPa)\n',PinGPa);
	fprintf(fileID,'Time for pressure calculation      :  %.3f (sec)\n',toc(t1));
		
	fclose(fileID);
end 

% Updating relaxation iteration number
S.S_Debug.relax(S.Relax_iter).relax_time = toc(tic_relax);
if (S.RelaxFlag || S.MDFlag)
	fprintf('\n Relaxation step number %d completed in %f s.\n',S.Relax_iter, S.S_Debug.relax(S.Relax_iter).relax_time);

	outfname = S.outfname;
	fileID = fopen(outfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',outfname);
	end
	fprintf(fileID,'Relax time                         :  %.3f (sec)\n',toc(tic_relax));
	fclose(fileID);

end
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

atompos_next = 0;
AtmForce = 0;

% Reshape the atomic positions for Relax/MD
if (S.RelaxFlag == 1 || S.MDFlag)
    if S.cell_typ <= 2
	    %atompos_next = transpose(S.lat_uvec) * S.Atoms';
	    atompos_next = coordinateTransformation(S, S.Atoms, 'noncart2cart_dis');
	    % update mapped Atomic position history, for keeping track of un-mapped
	    % positions that are used in charge extrapolation
	    S.Atoms_old = atompos_next;
	    atompos_next = reshape(atompos_next',[],1);
	    % -Forces
        AtmForce = reshape(S.force',[],1);
        AtmForce = - AtmForce;
	else
	    atompos_next = reshape(atom_pos',[],1);
        atompos_noncart = coordinateTransformation(S, atom_pos, 'cart2noncart_dis');
        %fmag_xy = sqrt(S.force(:,1).^2 + S.force(:,2).^2);
        %AtmForce = S.force;
        %AtmForce(:,1) = fmag_xy .* cos(atompos_noncart(:,2) + S.twist*atompos_noncart(:,3));
        %AtmForce(:,2) = fmag_xy .* sin(atompos_noncart(:,2) + S.twist*atompos_noncart(:,3));
        AtmForce = S.force';
		for i = 1:size(AtmForce,2)
			fac1 = (atompos_noncart(i,2) - S.Atoms(i,2))/S.L2;
			fac2 = (atompos_noncart(i,3) - S.Atoms(i,3))/S.L3;
			ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
			AtmForce(:,i) = ROT * AtmForce(:,i);
		end
		
		AtmForce = reshape(AtmForce,[],1);
        AtmForce = - AtmForce;
    end
end

S=phonon(S);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = initElectrondensity(S)
if ((S.elecgs_Count-S.StressCount+1) == 1)
	% use sum of the atomic charge densities as rho_guess
	S.rho = S.rho_at;
else
	% perform charge extrapolation
	if ((S.RelaxFlag == 1 || S.MDFlag) && (S.elecgs_Count-S.StressCount+1) >= 4)
		S.rho(:,1) = S.rho_at(:,1) + S.delta_rho_in_tp1;
		% update spin up/down densities
		if S.nspin ~= 1
			rho_mag = S.rho(:,2) - S.rho(:,3);
			S.rho(:,2) = (S.rho(:,1) + rho_mag) * 0.5;
			S.rho(:,3) = (S.rho(:,1) - rho_mag) * 0.5;
		end
	end
	% need to scale the density
	S.rho = S.rho * abs(S.NegCharge/dot(S.W,S.rho(:,1)));
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = check_atomlocation(S)
% map atom positions back to the domain if atoms are outside the domain in
% periodic directions or throw an error in Dirichlet directions
if S.BCx == 1
	if (sum(S.Atoms(:,1) >= S.L1 | S.Atoms(:,1) < 0) > 0)
		error('Atom coordinates in the first lattice vector direction are out of cell');
	end
else
	%S.Atoms(S.Atoms(:,1) >= S.L1,1) = S.Atoms(S.Atoms(:,1) >= S.L1,1) - S.L1;
	%S.Atoms(S.Atoms(:,1) < 0,1) = S.Atoms(S.Atoms(:,1) < 0,1) + S.L1;
	 S.Atoms(:,1) = mod(S.Atoms(:,1),S.L1);
end

if S.BCy == 1
	if (sum(S.Atoms(:,2) >= S.L2 | S.Atoms(:,2) < 0) > 0)
		error('Atom coordinates in the second lattice vector direction are out of cell');
	end
else
	%S.Atoms(S.Atoms(:,2) >= S.L2,2) = S.Atoms(S.Atoms(:,2) >= S.L2,2) - S.L2;
	%S.Atoms(S.Atoms(:,2) < 0,2) = S.Atoms(S.Atoms(:,2) < 0,2) + S.L2;
	S.Atoms(:,2) = mod(S.Atoms(:,2),S.L2);
end

if S.BCz == 1
	if (sum(S.Atoms(:,3) >= S.L3 | S.Atoms(:,3) < 0) > 0)
		error('Atom coordinates in the third lattice vector direction are out of cell');
	end
else
	%S.Atoms(S.Atoms(:,3) >= S.L3,3) = S.Atoms(S.Atoms(:,3) >= S.L3,3) - S.L3;
	%S.Atoms(S.Atoms(:,3) < 0,3) = S.Atoms(S.Atoms(:,3) < 0,3) + S.L3;
	S.Atoms(:,3) = mod(S.Atoms(:,3),S.L3);	
end

end










