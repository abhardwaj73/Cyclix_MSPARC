function [W] = intgWts_cychel(Nx,Ny,Nz,BCx,BCy,BCz,xin,S)
    W_x = ((xin + [0:Nx-1]*S.dx) * S.dx)';
    W_x(1) = 0.5*xin*S.dx + 0.125*S.dx*S.dx;
    W_x(Nx) = 0.5*(xin + (Nx-1)*S.dx)*S.dx - 0.125*S.dx*S.dx;

    W_y = ones(Ny,1)*S.dy;
    W_y(1) = W_y(1) * (1-BCy*0.5);
    W_y(Ny) = W_y(Ny) * (1-BCy*0.5);

    W_z = ones(Nz,1)*S.dz;
    W_z(1) = W_z(1) * (1-BCz*0.5);
    W_z(Nz) = W_z(Nz) * (1-BCz*0.5);

    W = kron(W_z,kron(W_y,W_x));
end