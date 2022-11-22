function [DOS, E] = eig2DOS(lambda_nk, N, sigma)
% EIG2DOS plots density of states (DOS) based on the eigenvalues
% DOE(e) = 2 \sum_{n=1}^{Ns} delta(lambda(n) -  e)
%
% @param lambda Eigenvalues
% @param n      Number of points in E
% @param sigma  Standard deviation
%
% @author       Qimen Xu <qimenxu@gatech.edu>

N_k = size(lambda_nk,1);

lambda = reshape(lambda_nk,[],1);

lambda_min = min(lambda);
lambda_max = max(lambda);
width = lambda_max - lambda_min;

%buf = max(5*sigma+0.1, width*0.1);
%buf = 5*sigma+width*0.01;
buf = 0;

Emin = lambda_min - buf;
Emax = lambda_max + buf;

%N = 5*length(lambda);

E = linspace(Emin, Emax, N);

DOS = sum(gauss_distribution(colminusrow(E, lambda), sigma), 2)/N_k;

end


function f = gauss_distribution(x, s)
% f(x,mu,s) = 1/(s*sqrt(2*pi)) * exp( -0.5 * ((x)/s)^2 )
p1 = -(x/s) .^ 2;
p2 = (s * sqrt(pi));
f = exp(p1)/p2; 
end


% tool function colminusrow
function xmy = colminusrow(x,y)
% A column vector x minus a row vector.
% In Matlab versions after R2018b, it's just x - y
[xx,yy] = ndgrid(x,y);
xmy = xx - yy;
end






