% Thu  6 Jan 15:53:31 CET 2022
% Karl Kästner, Berlin
%
%% normaliztation scale of the spatial bandpass density
%
%function [Sc] = spectral_density_bandpass_continuous_scale(fc,p)
function [Sc] = spectral_density_bandpass_continuous_scale(fc,p,numeric,q)
	if (nargin()<3)
		numeric = false;
	end
	if (0) %issym(p) || (p<43 && ~numeric) )
		% TODO this is not any more correct
		% Sc = (4*gamma(4*p-1))./(16.^p.*fc.*gamma((4*p-1)/2).^2);
	else
		% avoid overflow
		% TODO these limits are for q = 1
		tol = 1e-5^(1/(4*p));
		fl = fc*(1-sqrt(1 - tol^2))/tol;
		fr = fc*(1+sqrt(1 - tol^2))/tol;
		Sc = 1./quad(@(fx) spectral_density_bandpass_continuous(fx,fc,p,-1,q),fl,fr);
	end
end

