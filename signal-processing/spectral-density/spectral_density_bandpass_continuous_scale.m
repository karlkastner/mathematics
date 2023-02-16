% Thu  6 Jan 15:53:31 CET 2022
% Karl Kästner, Berlin
%
%% normaliztation scale of the spatial bandpass density
%
% function [Sc] = spectral_density_bandpass_continuous_scale(fc,p)
function [Sc] = spectral_density_bandpass_continuous_scale(fc,p,pp,numeric)
	if (nargin()<3)
		pp = [];
	end
	if (nargin()<4)
		numeric = false;
	end
%	else
	if (isempty(pp) && ~numeric)
		kc = 2*pi*fc;
		% IS = kc./(2*pi) * 4.^(1-p)./(2*p-1).*pi.*p.*binom(2*p-1,p-1)
		IS = kc.*2.^(1-2*p).*binom(2*(p-1),p-1);
		if (isnan(IS))
			% stirlings approximation
			IS = kc./sqrt(4*pi*(p-1));
		end
		if (0 == IS)
			Sc = realmax;
		else
			Sc = 1./IS;
		end
	else
		% avoid overflow
		% TODO these limits are for q = 1
		tol = 1e-5^(1/(4*p));
		fl = fc*(1-sqrt(1 - tol^2))/tol;
		fr = fc*(1+sqrt(1 - tol^2))/tol;
		Sc = 1./quad(@(fx) spectral_density_bandpass_continuous(fx,fc,p,0,pp),fl,fr);
	end
	if (~issym(p))
		Sc(p<0.5) = NaN;
	end
end

