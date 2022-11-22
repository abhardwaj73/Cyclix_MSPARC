function dd = calculateDistance_cychel(X,Y,Z,X_ref,Y_ref,Z_ref,S);
    pos_node_ch = [X(:),Y(:),Z(:)];
    pos_node_cart = coordinateTransformation_cychel(S,pos_node_ch,'noncart2cart_dis');
    pos_atm_cart = coordinateTransformation_cychel(S,[X_ref,Y_ref,Z_ref],'noncart2cart_dis');
    dd = reshape(sum(bsxfun(@minus,pos_node_cart,pos_atm_cart).^2,2).^0.5,size(X,1),size(X,2),size(X,3));
end