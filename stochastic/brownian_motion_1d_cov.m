% 2022-06-13 11:33:33.337417529 +0200
%
% generate 1d noise with the density of (fractional) brownian motion
% H = 1/2 : brownian noise (wiener process)
%
% Kaplan 1996
% c.f. Danudirdjo 2011
% 
function [y,C] = brownian_motion_1d_cov(n,H)
	if (length(n)<2)
		n(2)=1;
	end
	if (nargin()<2)
		H = 1/2;
	end

	e = randn(n(1),n(2));
	y = zeros(n(1));
	id = (1:n);
	% covariance of the first row without circular boundary conditions
	%c = abs(id_+1).^h - 2*id_.^h + abs(id_-1).^h;
	%c = 1/2*(abs(k+1) - 2*k + abs(k-1));
	% c.f. Danudirdjo 2011
	if (H~=1/2)
		C = 1/2*(id.^(2*H) + id'.^(2*H) - abs(id-id').^(2*H));
	else
		C = 1/2*(id + id' - abs(id-id'));
	end
	y = sqrtm(C)*e;
end % brownian_motion_1d_cov

