function S = gradIndicesValues_cychel(S)
Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;
N = S.N;
n0 = S.FDn;
w1 = S.w1;
dx = S.dx;
dy = S.dy;
dz = S.dz;

% Initial number of non-zeros: including ghost nodes
nnz_count = 4*n0*N ;

% Row numbers and non-zero values
I = zeros(nnz_count,1) ;
V1 = zeros(nnz_count,1) ;
V2 = zeros(nnz_count,1) ;
V3 = zeros(nnz_count,1) ;

% Indices of the columns
II = zeros(nnz_count,1);
JJ = zeros(nnz_count,1);
KK = zeros(nnz_count,1);

% Gradient along cartesian x and y directions
row_count = 1;
count = 1 ;
for kk=1:Nz
    z = (kk-1)*dz;
    for jj=1:Ny
        y = (jj-1)*dy;
        for ii=1:Nx
            x = S.xin + (ii-1)*dx;
            pos_node_cart = coordinateTransformation(S,[x,y,z],'noncart2cart_dis');
            wt1_d1 = pos_node_cart(1)/(x*dx);
            wt1_d2 = -pos_node_cart(2)/(x*x*dy);
            wt2_d1 = pos_node_cart(2)/(x*dx);
            wt2_d2 = pos_node_cart(1)/(x*x*dy);
            % off-diagonal elements
            for p=1:n0
                % ii+p
                I(count) = row_count; II(count) = ii+p; JJ(count) = jj; KK(count) = kk;
                V1(count) = wt1_d1*w1(p+1);
                V2(count) = wt2_d1*w1(p+1);
                count = count + 1;
                % ii-p
                I(count) = row_count; II(count) = ii-p ; JJ(count) = jj; KK(count) = kk;
                V1(count) = -wt1_d1*w1(p+1);
                V2(count) = -wt2_d1*w1(p+1);
                count = count + 1;
                % jj+p
                I(count) = row_count; II(count) = ii; JJ(count) = jj+p; KK(count) = kk;
                V1(count) = wt1_d2*w1(p+1);
                V2(count) = wt2_d2*w1(p+1);
                count = count + 1;
                % jj-p
                I(count) = row_count; II(count) = ii; JJ(count) = jj-p; KK(count) = kk;
                V1(count) = -wt1_d2*w1(p+1);
                V2(count) = -wt2_d2*w1(p+1);
                count = count + 1;
            end
            row_count = row_count+1;
        end
    end
end

isIn = (II >= 1) & (II <= Nx);
I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V1 = V1(isIn); V2 = V2(isIn);

if S.BCy == 1
    % Removing outside domain entries (for periodic code this is unnecessary)
    isIn = (JJ >= 1) & (JJ <= Ny);
    I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V1 = V1(isIn); V2 = V2(isIn);
elseif S.BCy == 0
    S.G_JOutl_1 = (JJ<1); S.G_JOutr_1 = (JJ>Ny); % Warning: Assumed influence of only neighboring cells
    JJ = mod(JJ+(Ny-1),Ny)+1;
end

% Getting linear indices of the columns
S.G_J_1 = (KK-1)*Nx*Ny + (JJ-1)*Nx + II;
S.G_I_1 = I;
S.G_V_1 = V1;

S.G_J_2 = S.G_J_1;
S.G_I_2 = I;
S.G_V_2 = V2;

% Gradient along cartesian z directions
row_count = 1;
count = 1 ;
wt3_d2 = -S.twist/dy;
wt3_d3 = 1/dz;
for kk=1:Nz
    for jj=1:Ny
        for ii=1:Nx
            % off-diagonal elements
            for p=1:n0
                % jj+p
                I(count) = row_count; II(count) = ii; JJ(count) = jj+p; KK(count) = kk;
                V3(count) = wt3_d2*w1(p+1);
                count = count + 1;
                % jj-p
                I(count) = row_count; II(count) = ii; JJ(count) = jj-p; KK(count) = kk;
                V3(count) = -wt3_d2*w1(p+1);
                count = count + 1;
                % kk+p
                I(count) = row_count; II(count) = ii; JJ(count) = jj; KK(count) = kk+p;
                V3(count) = wt3_d3*w1(p+1);
                count = count + 1;
                % kk-p
                I(count) = row_count; II(count) = ii; JJ(count) = jj; KK(count) = kk-p;
                V3(count) = -wt3_d3*w1(p+1);
                count = count + 1;
                
            end
            row_count = row_count+1;
        end
    end
end

if S.BCy == 1
    % Removing outside domain entries (for periodic code this is unnecessary)
    isIn = (JJ >= 1) & (JJ <= Ny);
    I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V3 = V3(isIn);
elseif S.BCy == 0
    S.G_JOutl_2 = (JJ<1); S.G_JOutr_2 = (JJ>Ny); % Warning: Assumed influence of only neighboring cells
    JJ = mod(JJ+(Ny-1),Ny)+1;
end
if S.BCz == 1
    % Removing outside domain entries (for periodic code this is unnecessary)
    isIn = (KK >= 1) & (KK <= Nz);
    I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V3 = V3(isIn);
elseif S.BCz == 0
    S.G_JOutl_3 = (KK<1); S.G_JOutr_3 = (KK>Nz); % Warning: Assumed influence of only neighboring cells
    KK = mod(KK+(Nz-1),Nz)+1;
end

% Getting linear indices of the columns
S.G_J_3 = (KK-1)*Nx*Ny + (JJ-1)*Nx + II;
S.G_I_3 = I;
S.G_V_3 = V3;
end
