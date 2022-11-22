function [DX_x,DX_y,DX_z] = dpseudopot_cychel(X,II,JJ,KK,XX,YY,ZZ,S)
	DX_x = zeros(size(X));
    DX_y = zeros(size(X));
    DX_z = zeros(size(X));

    % Convert to cartesian coordinates
    pos_node_cart = coordinateTransformation_cychel(S,[XX(:),YY(:),ZZ(:)],'noncart2cart_dis');
    XX_cart = reshape(pos_node_cart(:,1),size(XX,1),size(XX,2),size(XX,3));
    YY_cart = reshape(pos_node_cart(:,2),size(YY,1),size(YY,2),size(YY,3));
    ZZ_cart = reshape(pos_node_cart(:,3),size(ZZ,1),size(ZZ,2),size(ZZ,3));
    
    % Coefficients in the gradient
    c11 = XX_cart./(XX*S.dx);
    c12 = -YY_cart./(XX.*XX*S.dy);
    c21 = YY_cart./(XX*S.dx);
    c22 = XX_cart./(XX.*XX*S.dy);
    c31 = -S.twist/S.dy;
    c32 = 1/S.dz ;

    for p = 1:S.FDn
    	coef = S.w1(p+1);   
		DX_x(II,JJ,KK) = DX_x(II,JJ,KK) + c11*coef.*(X(II+p,JJ,KK)-X(II-p,JJ,KK)) ...
                                        + c12*coef.*(X(II,JJ+p,KK)-X(II,JJ-p,KK)) ;
		DX_y(II,JJ,KK) = DX_y(II,JJ,KK) + c21*coef.*(X(II+p,JJ,KK)-X(II-p,JJ,KK)) ...
                    					+ c22*coef.*(X(II,JJ+p,KK)-X(II,JJ-p,KK)) ;
		DX_z(II,JJ,KK) = DX_z(II,JJ,KK) + c31*coef.*(X(II,JJ+p,KK)-X(II,JJ-p,KK)) ...
	           							+ c32*coef.*(X(II,JJ,KK+p)-X(II,JJ,KK-p)) ;
    end

end
