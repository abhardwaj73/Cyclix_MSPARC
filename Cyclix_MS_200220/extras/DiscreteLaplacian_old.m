function DL = DiscreteLaplacian_old(S,kpt)

Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;
N  = S.N; n0 = S.FDn;
w1 = S.w1;
w2 = S.w2;
dx = S.dx;
dy = S.dy;
dz = S.dz;


% Phase factor
% PhaseFactor = exp(-pi*1i*nu); % it is equivalent to (-1)^nu

% Initial number of non-zeros: including ghost nodes
nnzCount = (3 * (2 * n0) + 1) * N;

% Row and column indices and the corresponding non-zero values
% used to generate sparse matrix DL s.t. DL(I(k),J(k)) = V(k)
I = zeros(nnzCount,1);
V = zeros(nnzCount,1);

% Indices of the columns J = II*Ntheta*(Nz+1) + JJ*(Nz+1) + KK + 1;
II = zeros(nnzCount,1);
JJ = zeros(nnzCount,1);
KK = zeros(nnzCount,1);

rowCount = 1;
count = 1;
dx2 = dx^2;
dy2 = dy^2;
dz2 = dz^2;

% Find non-zero entries that use forward difference
for kk = 1:Nz
    for jj = 1:Ny
        for ii = 1:Nx
            % diagonal element
            I(count) = rowCount; II(count) = ii; JJ(count) = jj; KK(count) = kk;
            V(count) = S.lapc_T(1,1)*w2(1)/(dx2) + S.lapc_T(2,2)*w2(1)/(dy2) + S.lapc_T(3,3)*w2(1)/(dz2);
            count = count + 1;
            % off-diagonal elements
            for q = 1:n0
                % ii + q
                I(count) = rowCount; II(count) = ii+q; JJ(count) = jj; KK(count) = kk;
                V(count) = S.lapc_T(1,1)*w2(q+1)/(dx2);
                count = count + 1;
                % ii - q
                I(count) = rowCount; II(count) = ii-q; JJ(count) = jj; KK(count) = kk;
                V(count) = S.lapc_T(1,1)*w2(q+1)/(dx2);
                count = count + 1;
                % jj + q
                I(count) = rowCount; II(count) = ii; JJ(count) = jj+q; KK(count) = kk;
                V(count) = S.lapc_T(2,2)*w2(q+1)/(dy2);
                count = count + 1;
                % jj - q
                I(count) = rowCount; II(count) = ii; JJ(count) = jj-q; KK(count) = kk;
                V(count) = S.lapc_T(2,2)*w2(q+1)/(dy2);
                count = count + 1;
                % kk + q
                I(count) = rowCount; II(count) = ii; JJ(count) = jj; KK(count) = kk+q;
                V(count) = S.lapc_T(3,3)*w2(q+1)/(dz2);
                count = count + 1;
                % kk - q
                I(count) = rowCount; II(count) = ii; JJ(count) = jj; KK(count) = kk-q;
                V(count) = S.lapc_T(3,3)*w2(q+1)/(dz2);
                count = count + 1;
            end
            rowCount = rowCount + 1;
        end
    end
end


if (S.BC == 1)
    isIn = (II >= 1) & (II <= Nx) & (JJ >= 1) & (JJ <= Ny) & (KK >= 1) & (KK <= Nz);
    I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V = V(isIn);
elseif (S.BC == 2)
    isIn11 = (II<1);isIn12 = (II>Nx);isIn21 = (JJ<1);isIn22 = (JJ>Ny);isIn31 = (KK<1);isIn32 = (KK>Nz); % Warning: Assumed influence of only neighboring cells
    V(isIn11) = V(isIn11)*exp(-1i*2*pi*kpt(1)); V(isIn12) = V(isIn12)*exp(1i*2*pi*kpt(1));
    V(isIn21) = V(isIn21)*exp(-1i*2*pi*kpt(2)); V(isIn22) = V(isIn22)*exp(1i*2*pi*kpt(2));
    V(isIn31) = V(isIn31)*exp(-1i*2*pi*kpt(3)); V(isIn32) = V(isIn32)*exp(1i*2*pi*kpt(3));
    II = mod(II-1+S.Nx,S.Nx)+1; JJ = mod(JJ-1+S.Ny,S.Ny)+1; KK = mod(KK-1+S.Nz,S.Nz)+1;
elseif (S.BC == 3)
    isIn = (KK >= 1) & (KK <= Nz);
    I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V = V(isIn);
    isIn11 = (II<1);isIn12 = (II>Nx);isIn21 = (JJ<1);isIn22 = (JJ>Ny);
    V(isIn11) = V(isIn11)*exp(-1i*2*pi*kpt(1)); V(isIn12) = V(isIn12)*exp(1i*2*pi*kpt(1));
    V(isIn21) = V(isIn21)*exp(-1i*2*pi*kpt(2)); V(isIn22) = V(isIn22)*exp(1i*2*pi*kpt(2));
    II = mod(II-1+S.Nx,S.Nx)+1; JJ = mod(JJ-1+S.Ny,S.Ny)+1;
elseif (S.BC == 4)
    isIn = (II >= 1) & (II <= Nx) & (JJ >= 1) & (JJ <= Ny);
    I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V = V(isIn);
    isIn31 = (KK<1);isIn32 = (KK>Nz);
    V(isIn31) = V(isIn31)*exp(-1i*2*pi*kpt(3)); V(isIn32) = V(isIn32)*exp(1i*2*pi*kpt(3));
    KK = mod(KK-1+S.Nz,S.Nz)+1;
end


% Getting linear indices of the columns
J = (KK-1)*Nx*Ny + (JJ-1)*Nx + II;

% Create discretized Laplacian
DL = sparse(I,J,V,N,N);

%  Add Mixed Derivative Componenets
%-----------------------------------

dxdy = S.dx*S.dy;
dydz = S.dy*S.dz;
dzdx = S.dz*S.dx;
% Create 1D gradient in all directions
%---------------------------------------

if (S.cell_typ == 2)
    
    nnz_x = 2*n0*Nx;
    G = zeros(nnz_x,1);
    R = zeros(nnz_x,1);
    A = zeros(nnz_x,1);
    rowCount = 1;
    count = 1;
    
    % x-direction
    %-------------
    
    for ii = 1:Nx
        for q = 1:n0
            % ii + q
            G(count) = rowCount; R(count) = ii+q;
            A(count) = w1(q+1);
            count = count + 1;
            % ii - q
            G(count) = rowCount; R(count) = ii-q;
            A(count) = -w1(q+1);
            count = count + 1;
        end
        rowCount = rowCount + 1;
    end
    
    isIn11 = (R<1); isIn12 = (R>Nx);
    A(isIn11) = A(isIn11)*exp(-1i*2*pi*kpt(1)); A(isIn12) = A(isIn12)*exp(1i*2*pi*kpt(1));
    R = mod(R-1+Nx,Nx)+1;

    GRA_x = sparse(G,R,A,Nx,Nx);
    
    % y-direction
    %-------------
    
    nnz_y = 2*n0*Ny;
    G = zeros(nnz_y,1);
    R = zeros(nnz_y,1);
    A = zeros(nnz_y,1);
    count =1;
    rowCount =1;
    
    for jj = 1:Ny
        for q = 1:n0
            % jj + q
            G(count) = rowCount; R(count) = jj+q;
            A(count) = w1(q+1);
            count = count + 1;
            % jj - q
            G(count) = rowCount; R(count) = jj-q;
            A(count) = -w1(q+1);
            count = count + 1;
        end
        rowCount = rowCount + 1;
    end
    
    isIn21 = (R<1); isIn22 = (R>Ny);
     A(isIn21) = A(isIn21)*exp(-1i*2*pi*kpt(2)); A(isIn22) = A(isIn22)*exp(1i*2*pi*kpt(2));
     R = mod(R-1+Ny,Ny)+1;

    GRA_y = sparse(G,R,A,Ny,Ny);
    
    if(S.BC == 2)
        
        % z-direction
        %-------------
        
        nnz_z = 2*n0*Nz;
        G = zeros(nnz_z,1);
        R = zeros(nnz_z,1);
        A = zeros(nnz_z,1);
        count =1;
        rowCount =1;
        
        for kk = 1:Nz
            for q = 1:n0
                % kk + q
                G(count) = rowCount; R(count) = kk+q;
                A(count) = w1(q+1);
                count = count + 1;
                % kk - q
                G(count) = rowCount; R(count) = kk-q;
                A(count) = -w1(q+1);
                count = count + 1;
            end
            rowCount = rowCount + 1;
        end
        
        isIn31 = (R<1); isIn32 = (R>Nz);
        A(isIn31) = A(isIn31)*exp(-1i*2*pi*kpt(3)); A(isIn32) = A(isIn32)*exp(1i*2*pi*kpt(3));
        R = mod(R-1+Nz,Nz)+1;

        GRA_z = sparse(G,R,A,Nz,Nz);
    end
end

% Create matrices for mixed derivatives
if (S.cell_typ == 2) && (S.BC == 2)
    MDL = S.lapc_T(1,2)/dxdy * kron(speye(Nz),kron(GRA_y,GRA_x))  +  S.lapc_T(2,3)/dydz * kron(GRA_z,kron(GRA_y,speye(Nx))) + ...
          S.lapc_T(1,3)/dzdx * kron(GRA_z,kron(speye(Ny),GRA_x)) ;
    DL = DL + MDL;
elseif (S.cell_typ == 2) && (S.BC == 3)
     MDL = S.lapc_T(1,2)/dxdy * kron(speye(Nz),kron(GRA_y,GRA_x));
     DL = DL + MDL;
end
%   MDL = (S.lapc_T(1,2) + S.lapc_T(2,1))/dxdy * kron(speye(Nz),kron(GRA_y,GRA_x));
%   DL = DL + MDL;   