function b = pseudochargeDensity_atom_cychel(V,II,JJ,KK,xin,S)
    b = zeros(size(V));
    % Coefficients in the Laplacian
    c1 = 1/(S.dx*S.dx) ;
    c2 = 1/(S.dy*S.dy) ;
    c3 = 1/(S.dz*S.dz) ;
    c4 = (S.twist*S.twist)/(S.dy*S.dy);
    c5 = -(2*S.twist)/(S.dy*S.dz) ;
    c6 = 1/S.dx ;

    rr = xin + (II-1)*S.dx;
    rr = rr(:);
    coeffbb = S.w2(1)*(c1+c2./(rr.*rr)+c3+c4);
    b(II,JJ,KK) = bsxfun(@times,coeffbb,V(II,JJ,KK)) ;

    for p=1:S.FDn
        coeffLap    = S.w2(p+1);
        coeffgrad   = S.w1(p+1);
        b(II,JJ,KK) = b(II,JJ,KK)   + bsxfun(@times, ( c1                * coeffLap + (c6*coeffgrad)./rr ) , V(II+p,JJ,KK) ) ...
                                    + bsxfun(@times, ( c1                * coeffLap - (c6*coeffgrad)./rr ) , V(II-p,JJ,KK) ) ...
                                    + bsxfun(@times, (c2./(rr.*rr) + c4) * coeffLap                        , V(II,JJ+p,KK) ) ...
                                    + bsxfun(@times, (c2./(rr.*rr) + c4) * coeffLap                        , V(II,JJ-p,KK) ) ...
                                    +                  c3                * coeffLap                        * V(II,JJ,KK+p  )   ...
                                    +                  c3                * coeffLap                        * V(II,JJ,KK-p  );
                                    
        for q=1:S.FDn
            coeff=c5*coeffgrad*S.w1(q+1);
            b(II,JJ,KK) = b(II,JJ,KK) + coeff * V(II,JJ+p,KK+q) ...
                                      - coeff * V(II,JJ+p,KK-q) ...
                                      - coeff * V(II,JJ-p,KK+q) ...
                                      + coeff * V(II,JJ-p,KK-q) ;
        end
    end
end