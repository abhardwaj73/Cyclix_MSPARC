function DLX = lapVec_cychel(DL11,DL22,DL33,DG1,DG2,DG3,X,S)
% @ brief    Calculates laplacian vector product using 
%            Kronecker product method
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param DLii        Discrete laplacian component in 1D along ith direction
% @param DGi         Discrete gradient component in 1D along ith direction
% @param X           N x Ns matrix of eigenfunctions of Hamiltonian
% @param DLX         Discrete laplacian times X
% @2016-2019 (c) Georgia Institute of Technology.
% @ references
%              "On real-space Density Functional Theory for 
%               non-orthogonal crystal systems: Kronecker product 
%               formulation of the kinetic energy operator (Sharma et. al. 2018)"
%=====================================================================================
Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;
if S.cell_typ == 3
	X1 = reshape(X,Nx,Ny,[]);
	l_zs = size(X1,3);
	Hlx1 = zeros(Nx,Ny,l_zs);
	for i=1:l_zs
	    Hlx1(:,:,i) = (DL11+DG1)*X1(:,:,i) + S.R2inv*X1(:,:,i)*DL22.';
	end
	Hlx1 = reshape(Hlx1,Nx*Ny,Nz,[]);
	X2 = reshape(X,Nx*Ny,Nz,[]);
	l_s = size(X2,3);
	for i=1:l_s
	    Hlx1(:,:,i) = Hlx1(:,:,i) + X2(:,:,i)*DL33.';
	end

	DLX = reshape(Hlx1,[],l_s);
elseif (S.cell_typ == 4 || S.cell_typ == 5)
	X1 = reshape(X,Nx,Ny,[]);
	l_zs = size(X1,3);
	Hlx1 = zeros(Nx,Ny,l_zs);
	Hlx2 = zeros(Nx,Ny,l_zs);
	for i=1:l_zs
	    Hlx1(:,:,i) = (DL11+DG1)*X1(:,:,i) + X1(:,:,i)*DL22.' + S.R2inv*X1(:,:,i)*DL22.';
	    Hlx2(:,:,i) = X1(:,:,i)*DG2.';
	end
	Hlx1 = reshape(Hlx1,Nx*Ny,Nz,[]) ;
	Hlx2 = reshape(Hlx2,Nx*Ny,Nz,[]) ;

	X2 = reshape(X,Nx*Ny,Nz,[]);
	l_s = size(X2,3);
	for i=1:l_s
	    Hlx2(:,:,i) = X2(:,:,i)*DL33.' + Hlx2(:,:,i)*DG3.';
	end

	DLX = reshape((Hlx1 + Hlx2),[],l_s);
end