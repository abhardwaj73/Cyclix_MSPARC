function phi = poisson_RHS(rho,b,L1,L2,L3,Nx,Ny,Nz,Rin)
	% rho - electron density
	% b - total pseudocharge density
	% L1,L2,L3 - lengths of the domain as given in the .inpt file
	% Nx,Ny,Nz - fd nodes in three directions as given in the .inpt file
	% Rin - internal radius of the tube
	
	f = -4 * pi * (rho + b);
	
	%reshape rho to 3D 
	rho = reshape(f, Nx, Ny*Nz); % note here after rho = rho + b
	
	dx = L1 / Nx;
	dy = L2 / Ny;
	dz = L3 / Nz;

	%sum over Z direction, \int (\rho) dz / LZ
	r = [Rin + (0:Nx-1) * dx]'; 
	
	phi = (1/2/pi/L2/L3)*sum(sum(rho,2).*log(r)) * dx *dy * dz;
end