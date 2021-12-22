% Sat  3 Jul 12:17:56 CEST 2021
% mu : mean
% sd : relative standard deviation per unit meter
% dx : grid distance, over which is averaged
% dt : time step (only for variation in time, put 1 for perturbed coefficients)
%
% E[r]        = mu
% std(r)      = mu*sd/sqrt(dt*dx)
% skewness(r) = 2*s/sqrt(k)
% kurtosis(r) = 3 + 6*s^2/(dx*dy)
function r = spatialrnd(mu,sdrel1,dt,dx,dy,n1,n2)
	if (nargin()<6)
		n2 = 1;
	end
	k = dx*dy*dt/(sdrel1*sdrel1);
	r = mu*gamrnd(k,1/k,n1,n2);
end

