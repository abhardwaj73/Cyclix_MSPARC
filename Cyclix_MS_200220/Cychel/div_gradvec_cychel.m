function [M1,M2] = div_gradvec_cychel(S,kptvec)
% @ brief   Calculates equivalent gradient matrices that need to operate on F1 and F2 of F vector
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param Mi            matrices that will operate on F1 and F2 of F vector 
% @2016-2019 (c) Georgia Institute of Technology.
%===================================================================================


% F = [F1 F2 F3] is a vector
% dF1/dx1 + dF2/dx2 need to be calculated (cartessian gradients)

I = S.G_I_1;
J = S.G_J_1;

V1_1 = S.G_V_1; % coeeficients for d/dx1
V2_2 = S.G_V_2; % coefficients for d/dx2

if S.BCy == 0

    % multiply rotation factors to the fd coefficients for the nodes on the left of domain
    Outdoml = S.G_JOutl_1;
    V1_2l = S.G_V_1(Outdoml);
    V2_1l = S.G_V_2(Outdoml);

    L = S.L2;
    kpt = kptvec(2);
    if (kpt == 0)
        phasefac_l = 1.0;
        phasefac_r = 1.0;
    else
        phasefac_l = exp(-1i*kpt*L);
        phasefac_r = exp(1i*kpt*L);
    end

    % F1(R^-1 * x) = cos(theta) * F1(x) + sin(theta) * F2(x)
    V1_1(Outdoml) = V1_1(Outdoml) * cos(S.L2) * phasefac_l; % contribution from F1 to dF1/dx1
    V1_2l = V1_2l * sin(S.L2) * phasefac_l; % contribution from F2 to dF1/dx1

    % F2(R^-1 * x) = -sin(theta) * F1(x) + cos(theta) * F2(x)
    V2_1l = V2_1l * -sin(S.L2) * phasefac_l; % contribution from F1 to dF2/dx2
    V2_2(Outdoml) = V2_2(Outdoml) * cos(S.L2) * phasefac_l; % contribution from F2 to dF2/dx2
    
    % multiply rotation factors to the fd coefficients for the nodes on the right of domain
    Outdomr = S.G_JOutr_1;
    V1_2r = S.G_V_1(Outdomr);
    V2_1r = S.G_V_2(Outdomr);

    % F1(R * x) = cos(theta) * F1(x) - sin(theta) * F2(x)
    V1_1(Outdomr) = V1_1(Outdomr) * cos(S.L2) * phasefac_r; % contribution from F1 to dF1/dx1
    V1_2r = V1_2r * -sin(S.L2) * phasefac_r; % contribution from F2 to dF1/dx1

    % F2(R * x) = sin(theta) * F1(x) + cos(theta) * F2(x)
    V2_1r = V2_1r * sin(S.L2) * phasefac_r; % contribution from F1 to dF2/dx2
    V2_2(Outdomr) = V2_2(Outdomr) * cos(S.L2) * phasefac_r; % contribution from F2 to dF2/dx2

    V1_1(Outdoml) = V1_1(Outdoml) + V2_1l;
    V1_1(Outdomr) = V1_1(Outdomr) + V2_1r;
    V2_2(Outdoml) = V2_2(Outdoml) + V1_2l;
    V2_2(Outdomr) = V2_2(Outdomr) + V1_2r;
end

% Generate matrix that need to be multiplied by F1
M1 = sparse(I,J,V1_1,S.N,S.N);

% Generate matrix that need to be multiplied by F2
M2 = sparse(I,J,V2_2,S.N,S.N);

% if (div_dir == 1)
%     I = S.G_I_1;
%     J = S.G_J_1;
%     V1 = S.G_V_1;
%     V2 = S.G_V_2;
%     if S.BCy == 0
%         isOutl = S.G_JOutl_1;
%         isOutr = S.G_JOutr_1;
        
%         % boundary contribution from df/dx to d/dx(df/dx)
%         V1(isOutl) = V1(isOutl) * cos(S.L2);
%         V1(isOutr) = V1(isOutr) * cos(S.L2);

%         % boundary contribution from df/dx to d/dy(df/dy)
%         V2(isOutl) = V2(isOutl) * -sin(S.L2);      
%         V2(isOutr) = V2(isOutr) * sin(S.L2);
        
%         % collect components which are going to operate on df/dx
%         V1(isOutl) = V1(isOutl) + V2(isOutl);
%         V1(isOutr) = V1(isOutr) + V2(isOutr);
%     end
%     V = V1;
% elseif (div_dir == 2)
%     I = S.G_I_1;
%     J = S.G_J_1;
%     V1 = S.G_V_1;
%     V2 = S.G_V_2;
%     if S.BCy == 0
%         isOutl = S.G_JOutl_1;
%         isOutr = S.G_JOutr_1;

%         % boundary contribution from df/dy to d/dx(df/dx)
%         V1(isOutl) = V1(isOutl) * sin(S.L2);
%         V1(isOutr) = V1(isOutr) * -sin(S.L2);
        
%         % boundary contribution from df/dy to d/dy(df/dy)
%         V2(isOutl) = V2(isOutl) * cos(S.L2);
%         V2(isOutr) = V2(isOutr) * cos(S.L2);
        
%         % collect components which are going to operate on df/dy
%         V2(isOutl) = V2(isOutl) + V1(isOutl);
%         V2(isOutr) = V2(isOutr) + V1(isOutr);
%     end
%     V = V2;
% end

% DD = sparse(I,J,V,S.N,S.N);