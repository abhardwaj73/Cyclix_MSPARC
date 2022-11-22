function Y = coordinateTransformation_cychel(S,X,transfrm_typ)

if (strcmp(transfrm_typ,'noncart2cart_dis'))
    Y = zeros(size(X));
    Y(:,1) = X(:,1) .* cos(X(:,2) + S.twist * X(:,3) );
    Y(:,2) = X(:,1) .* sin(X(:,2) + S.twist * X(:,3) );
    Y(:,3) = X(:,3);
elseif (strcmp(transfrm_typ,'cart2noncart_dis'))
	Y = zeros(size(X));
	Y(:,1) = sqrt(X(:,1).^2 + X(:,2).^2);
	Y(:,2) = atan2(X(:,2),X(:,1)) - S.twist*X(:,3);
	Y(:,3) = X(:,3);
end