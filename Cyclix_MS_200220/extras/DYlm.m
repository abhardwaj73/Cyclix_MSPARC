function [DYlm_dx,DYlm_dy,DYlm_dz] = DYlm(lat_vec,X,Y,Z,l,m)
% function to perform the derivative of spherical harmonics

X_temp = lat_vec(1,1)*X + lat_vec(2,1)*Y + lat_vec(3,1)*Z;
Y_temp = lat_vec(1,2)*X + lat_vec(2,2)*Y + lat_vec(3,2)*Z;
Z_temp = lat_vec(1,3)*X + lat_vec(2,3)*Y + lat_vec(3,3)*Z;

X = X_temp;
Y = Y_temp;
Z = Z_temp;

r = sqrt(X.^2 + Y.^2 + Z.^2) ;

if (l == 0)
    % l=0
    C00 = 0.282094791773878; 	% 0.5*sqrt(1/pi)
    Ylm = C00 * ones(size(X));
    DYlm_dx = 0.0;
    DYlm_dy = 0.0;
    DYlm_dz = 0.0;
elseif (l == 1)
    % l=1
    C1m1 = 0.488602511902920; 	% sqrt(3/(4*pi))
    C10 = 0.488602511902920; 	% sqrt(3/(4*pi))
    C1p1 = 0.488602511902920; 	% sqrt(3/(4*pi))
    if(m==-1)
        Ylm = C1m1 * (Y ./ r);
        DYlm_dx = - C1m1 * (X.*Y)./(r.^3);
        DYlm_dy =   C1m1 *(r.*r - Y.*Y)./(r.^3);
        DYlm_dz = - C1m1 * (Z.*Y)./(r.^3);
    elseif(m==0)
        Ylm = C10 * (Z ./ r);
        DYlm_dx = - C10 * (X.*Z)./(r.^3);
        DYlm_dy = - C10 * (Z.*Y)./(r.^3);
        DYlm_dz =   C10 * (r.*r - Z.*Z)./(r.^3);
    elseif(m==1)
        Ylm = C1p1 * (X ./ r);
        DYlm_dx =   C1p1 * (r.*r - X.*X)./(r.^3);
        DYlm_dy = - C1p1 * (X.*Y)./(r.^3);
        DYlm_dz = - C1p1 * (Z.*X)./(r.^3);
    end
elseif(l == 2)
    % l=2
    C2m2 = 1.092548430592079; 	% 0.5*sqrt(15/pi)
    C2m1 = 1.092548430592079; 	% 0.5*sqrt(15/pi)
    C20 =  0.315391565252520; 	% 0.25*sqrt(5/pi)
    C2p1 = 1.092548430592079; 	% 0.5*sqrt(15/pi)
    C2p2 =  0.546274215296040;	% 0.25*sqrt(15/pi)
    if(m==-2)
        Ylm = C2m2*(X.*Y)./(r.*r);
        DYlm_dx = C2m2 * (-2*X.*X.*Y + r.*r.*Y)./(r.^4);
        DYlm_dy = C2m2 * (-2*X.*Y.*Y + r.*r.*X)./(r.^4);
        DYlm_dz = C2m2 * (-2*X.*Y.*Z)./(r.^4);
    elseif(m==-1)
        Ylm = C2m1*(Y.*Z)./(r.*r);
        DYlm_dx = C2m1 * (-2*X.*Y.*Z)./(r.^4);
        DYlm_dy = C2m1 * (-2*Z.*Y.*Y + r.*r.*Z)./(r.^4);
        DYlm_dz = C2m1 * (-2*Y.*Z.*Z + r.*r.*Y)./(r.^4);
    elseif(m==0)
        Ylm = C20*(-X.*X - Y.*Y + 2*Z.*Z)./(r.*r);
        DYlm_dx = C20 * (-2*X.*r.*r - 2*X.*(-X.*X - Y.*Y + 2*Z.*Z))./(r.^4);
        DYlm_dy = C20 * (-2*Y.*r.*r - 2*Y.*(-X.*X - Y.*Y + 2*Z.*Z))./(r.^4);
        DYlm_dz = C20 * (4*Z.*r.*r - 2*Z.*(-X.*X - Y.*Y + 2*Z.*Z))./(r.^4);
    elseif(m==1)
        Ylm = C2p1*(Z.*X)./(r.*r);
        DYlm_dx = C2p1 * (-2*Z.*X.*X + r.*r.*Z)./(r.^4);
        DYlm_dy = C2p1 * (-2*X.*Y.*Z)./(r.^4);
        DYlm_dz = C2p1 * (-2*X.*Z.*Z + r.*r.*X)./(r.^4);
    elseif(m==2)
        Ylm = C2p2*(X.*X - Y.*Y)./(r.*r);
        DYlm_dx = C2p2 * (2*X.*r.*r - 2*X.*(X.*X - Y.*Y))./(r.^4);
        DYlm_dy = C2p2 * (-2*Y.*r.*r - 2*Y.*(X.*X - Y.*Y))./(r.^4);
        DYlm_dz = C2p2 * (-2*Z.*(X.*X - Y.*Y))./(r.^4);
        
    end
elseif(l == 3)
    % l=3
    C3m3 =  0.590043589926644;	% 0.25*sqrt(35/(2*pi))
    C3m2 = 2.890611442640554;	% 0.5*sqrt(105/(pi))
    C3m1 = 0.457045799464466;	% 0.25*sqrt(21/(2*pi))
    C30 =  0.373176332590115; 	% 0.25*sqrt(7/pi)
    C3p1 =  0.457045799464466; 	% 0.25*sqrt(21/(2*pi))
    C3p2 = 1.445305721320277; 	% 0.25*sqrt(105/(pi))
    C3p3 = 0.590043589926644; 	% 0.25*sqrt(35/(2*pi))
    if(m==-3)
        Ylm = C3m3*(3*X.*X - Y.*Y).*Y./(r.*r.*r);
        DYlm_dx = C3m3*(6*X.*Y.*r.*r - 3*X.*Y.*(3*X.*X - Y.*Y))./(r.^5);
        DYlm_dy = C3m3*((3*X.*X - 3*Y.*Y).*r.*r - 3*Y.*Y.*(3*X.*X - Y.*Y))./(r.^5);
        DYlm_dz = C3m3*(-3*Z.*Y.*(3*X.*X - Y.*Y))./(r.^5);
    elseif(m==-2)
        Ylm = C3m2*(X.*Y.*Z)./(r.*r.*r);
        DYlm_dx = C3m2*(Y.*Z.*r.*r - 3*X.*X.*Y.*Z)./(r.^5);
        DYlm_dy = C3m2*(X.*Z.*r.*r - 3*X.*Y.*Y.*Z)./(r.^5);
        DYlm_dz = C3m2*(X.*Y.*r.*r - 3*X.*Y.*Z.*Z)./(r.^5);
    elseif(m==-1)
        Ylm = C3m1*Y.*(4*Z.*Z - X.*X - Y.*Y)./(r.*r.*r);
        DYlm_dx = C3m1*(-2*X.*Y.*r.*r - 3*X.*Y.*(-X.*X - Y.*Y + 4*Z.*Z))./(r.^5);
        DYlm_dy = C3m1*(-(X.*X + 3*Y.*Y -4*Z.*Z).*r.*r - 3*Y.*Y.*(-X.*X - Y.*Y + 4*Z.*Z))./(r.^5);
        DYlm_dz = C3m1*(8*Z.*Y.*r.*r - 3*Y.*Z.*(-X.*X - Y.*Y + 4*Z.*Z))./(r.^5);
    elseif(m==0)
        Ylm = C30*Z.*(2*Z.*Z - 3*X.*X - 3*Y.*Y)./(r.*r.*r);
        DYlm_dx = C30*(-6*X.*Z.*r.*r - 3*X.*Z.*(-3*X.*X - 3*Y.*Y + 2*Z.*Z))./(r.^5);
        DYlm_dy = C30*(-6*Y.*Z.*r.*r - 3*Y.*Z.*(-3*X.*X - 3*Y.*Y + 2*Z.*Z))./(r.^5);
        DYlm_dz = C30*(-(3*X.*X + 3*Y.*Y - 6*Z.*Z).*r.*r - 3*Z.*Z.*(-3*X.*X - 3*Y.*Y + 2*Z.*Z))./(r.^5);
    elseif(m==1)
        Ylm = C3p1*X.*(4*Z.*Z - X.*X - Y.*Y)./(r.*r.*r);
        DYlm_dx = C3p1*(-(3*X.*X + Y.*Y -4*Z.*Z).*r.*r - 3*X.*X.*(-X.*X - Y.*Y + 4*Z.*Z))./(r.^5);
        DYlm_dy = C3p1*(-2*X.*Y.*r.*r - 3*X.*Y.*(-X.*X - Y.*Y + 4*Z.*Z))./(r.^5);
        DYlm_dz = C3p1*(8*Z.*X.*r.*r - 3*X.*Z.*(-X.*X - Y.*Y + 4*Z.*Z))./(r.^5);
    elseif(m==2)
        Ylm = C3p2*Z.*(X.*X - Y.*Y)./(r.*r.*r);
        DYlm_dx = C3p2*(2*X.*Z.*r.*r - 3*X.*Z.*(X.*X - Y.*Y))./(r.^5);
        DYlm_dy = C3p2*(-2*Z.*Y.*r.*r - 3*Y.*Z.*(X.*X - Y.*Y))./(r.^5);
        DYlm_dz = C3p2*((X.*X - Y.*Y).*r.*r - 3*Z.*Z.*(X.*X - Y.*Y))./(r.^5);
    elseif(m==3)
        Ylm = C3p3*X.*(X.*X-3*Y.*Y)./(r.*r.*r);
        DYlm_dx = C3p3*((3*X.*X - 3*Y.*Y).*r.*r - 3*X.*X.*(X.*X - 3*Y.*Y))./(r.^5);
        DYlm_dy = C3p3*(-6*X.*Y.*r.*r - 3*X.*Y.*(X.*X - 3*Y.*Y))./(r.^5);
        DYlm_dz = C3p3*(-3*Z.*X.*(X.*X - 3*Y.*Y))./(r.^5);
    end
elseif(l == 4)
    if(m==-4)
        Ylm=(3.0/4.0)*sqrt(35.0/pi)*(X.*Y.*(X.*X-Y.*Y))./(r.^4);
    elseif(m==-3)
        Ylm=(3.0/4.0)*sqrt(35.0/(2.0*pi))*(3.0*X.*X-Y.*Y).*Y.*Z./(r.^4);
    elseif(m==-2)
        Ylm=(3.0/4.0)*sqrt(5.0/pi)*X.*Y.*(7.0*Z.*Z-r.*r)./(r.^4);
    elseif(m==-1)
        Ylm=(3.0/4.0)*sqrt(5.0/(2.0*pi))*Y.*Z.*(7.0*Z.*Z-3.0*r.*r)./(r.^4);
    elseif(m==0)
        Ylm=(3.0/16.0)*sqrt(1.0/pi)*(35.0*Z.^4-30.0*Z.*Z.*r.*r+3.0*r.^4)./(r.^4);
    elseif(m==1)
        Ylm=(3.0/4.0)*sqrt(5.0/(2.0*pi))*X.*Z.*(7.0*Z.*Z-3.0*r.*r)./(r.^4);
    elseif(m==2)
        Ylm=(3.0/8.0)*sqrt(5.0/(pi))*(X.*X-Y.*Y).*(7.0*Z.*Z-r.*r)./(r.^4);
    elseif(m==3)
        Ylm=(3.0/4.0)*sqrt(35.0/(2.0*pi))*(X.*X-3.0*Y.*Y).*X.*Z./(r.^4);
    elseif(m==4)
        Ylm=(3.0/16.0)*sqrt(35.0/pi)*(X.*X.*(X.*X-3.0*Y.*Y) - Y.*Y.*(3.0*X.*X-Y.*Y))./(r.^4);
    end
elseif(l == 5)
    p = sqrt(X.*X+Y.*Y);
    if(m==-5)
        Ylm = (3.0*sqrt(2*77/pi)/32.0)*(8.0*X.^4.*Y-4.0*X.*X.*Y.^3 + 4.0*Y.^5-3.0*Y.*p.^4)./(r.^5);
    elseif(m==-4)
        Ylm = (3.0/16.0)*sqrt(385.0/pi)*(4.0*X.^3.*Y - 4.0*X.*Y.^3).*Z./(r.^5);
    elseif(m==-3)
        Ylm = (sqrt(2.0*385.0/pi)/32.0)*(3.0*Y.*p.*p - 4.0*Y.^3).*(9*Z.*Z-r.*r)./(r.^5);
    elseif(m==-2)
        Ylm = (1.0/8.0)*sqrt(1155.0/pi)*2.0*X.*Y.*(3.0*Z.^3-Z.*r.*r)./(r.^5);
    elseif(m==-1)
        Ylm = (1.0/16.0)*sqrt(165.0/pi)*Y.*(21.0*Z.^4 - 14.0*r.*r.*Z.*Z+r.^4)./(r.^5);
    elseif(m==0)
        Ylm = (1.0/16.0)*sqrt(11.0/pi)*(63.0*Z.^5 -70.0*Z.^3.*r.*r + 15.0*Z.*r.^4)./(r.^5);
    elseif(m==1)
        Ylm = (1.0/16.0)*sqrt(165.0/pi)*X.*(21.0*Z.^4 - 14.0*r.*r.*Z.*Z+r.^4)./(r.^5);
    elseif(m==2)
        Ylm = (1.0/8.0)*sqrt(1155.0/pi)*(X.*X-Y.*Y).*(3.0*Z.^3 - r.*r.*Z)./(r.^5);
    elseif(m==3)
        Ylm = (sqrt(2.0*385.0/pi)/32.0)*(4.0*X.^3-3.0*p.*p.*X).*(9.0*Z.*Z-r.*r)./(r.^5);
    elseif(m==4)
        Ylm = (3.0/16.0)*sqrt(385.0/pi)*(4.0*(X.^4+Y.^4)-3.0*p.^4).*Z./(r.^5);
        %Ylm = (3.0/16.0)*sqrt(385.0/pi)*(X.^4+Y.^4-6*X.^2.*Y.^2).*Z./(r.^5);
    elseif(m==5)
        Ylm = (3.0*sqrt(2.0)/32.0)*sqrt(77.0/pi)*(4.0*X.^5 + 8.0*X.*Y.^4 -4.0*X.^3.*Y.*Y -3.0*X.*p.^4)./(r.^5);
    end
elseif(l == 6)
    p = sqrt(X.*X+Y.*Y);
    if(m==-6)
        Ylm = (sqrt(2.0*3003.0/pi)/64.0)*(12.0*X.^5.*Y+12.0*X.*Y.^5 - 8.0*X.^3.*Y.^3-6.0*X.*Y.*p.^4)./(r.^6);
    elseif(m==-5)
        Ylm = (3.0/32.0)*sqrt(2.0*1001.0/pi)*(8.0*X.^4.*Y - 4.0*X.*X.*Y.^3 + 4.0*Y.^5 -3.0*Y.*p.^4).*Z./(r.^6);
    elseif(m==-4)
        Ylm = (3.0/32.0)*sqrt(91.0/pi)*(4.0*X.^3.*Y -4.0*X.*Y.^3).*(11.0*Z.*Z-r.*r)./(r.^6);
    elseif(m==-3)
        Ylm = (sqrt(2.0*1365.0/pi)/32.0)*(-4.0*Y.^3 + 3.0*Y.*p.*p).*(11.0*Z.^3 - 3.0*Z.*r.*r)./(r.^6);
    elseif(m==-2)
        Ylm = (sqrt(2.0*1365/pi)/64.0)*(2.0*X.*Y).*(33.0*Z.^4-18.0*Z.*Z.*r.*r + r.^4)./(r.^6);
    elseif(m==-1)
        Ylm = (sqrt(273.0/pi)/16.0)*Y.*(33.0*Z.^5-30.0*Z.^3.*r.*r +5.0*Z.*r.^4)./(r.^6);
    elseif(m==0)
        Ylm = (sqrt(13.0/pi)/32.0)*(231.0*Z.^6-315*Z.^4.*r.*r + 105.0*Z.*Z.*r.^4 -5.0*r.^6)./(r.^6);
    elseif(m==1)
        Ylm = (sqrt(273.0/pi)/16.0)*X.*(33.0*Z.^5-30.0*Z.^3.*r.*r +5*Z.*r.^4)./(r.^6);
    elseif(m==2)
        Ylm = (sqrt(2.0*1365/pi)/64.0)*(X.*X-Y.*Y).*(33.0*Z.^4 - 18.0*Z.*Z.*r.*r + r.^4)./(r.^6);
    elseif(m==3)
        Ylm = (sqrt(2.0*1365.0/pi)/32.0)*(4.0*X.^3 -3.0*X.*p.*p).*(11.0*Z.^3 - 3.0*Z.*r.*r)./(r.^6);
    elseif(m==4)
        Ylm = (3.0/32.0)*sqrt(91.0/pi)*(4.0*X.^4+4.0*Y.^4 -3.0*p.^4).*(11.0*Z.*Z -r.*r)./(r.^6);
    elseif(m==5)
        Ylm = (3.0/32.0)*sqrt(2.0*1001.0/pi)*(4.0*X.^5 + 8.0*X.*Y.^4-4.0*X.^3.*Y.*Y-3.0*X.*p.^4).*Z./(r.^6);
    elseif(m==6)
        Ylm = (sqrt(2.0*3003.0/pi)/64.0)*(4.0*X.^6-4.0*Y.^6 +12.0*X.*X.*Y.^4-12.0*X.^4.*Y.*Y + 3.0*Y.*Y.*p.^4-3.0*X.*X.*p.^4)./(r.^6);
    end
end
if (l > 0)
    Ylm(abs(r) < 1e-10) = 0;
    DYlm_dx(abs(r) < 1e-10) = 0.0;
    DYlm_dy(abs(r) < 1e-10) = 0.0;
    DYlm_dz(abs(r) < 1e-10) = 0.0;
end

end
