% Produces unit lattice vectors from angles between edges
%----------------------------------------------------------
clc
clear

alpha = 80; % angle(in degree) between x-axis & z-axis
beta = 86.595132678838965; % angle between y-axis and z-axis
gamma = 70; % angles between x-axis and y-axis

uvec = zeros(3,3);
uvec(1,:) = [1 0 0];
uvec(2,:) = [cosd(gamma) sind(gamma) 0];
uvec(3,:) = [cosd(alpha) (cosd(beta) - cosd(alpha)*cosd(gamma))/sind(gamma) sqrt(sind(alpha)^2 - ((cosd(beta)-cosd(alpha)*cosd(gamma)) / sind(gamma))^2)];

fprintf('%.15f %.15f %.15f\n%.15f %.15f %.15f\n%.15f %.15f %.15f\n',uvec(1,1),uvec(1,2),uvec(1,3),uvec(2,1),uvec(2,2),uvec(2,3),uvec(3,1),uvec(3,2),uvec(3,3));
fprintf('norm : %.15f %.15f %.15f \n',norm(uvec(1,:)),norm(uvec(3,:)),norm(uvec(3,:)));


grad_T = inv(uvec') ;
lapc_T = grad_T * grad_T'