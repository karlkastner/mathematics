% Mon  6 Sep 10:10:13 CEST 2021
% function S = spreactral_density_lorentzian(fx,varargin)
function [S] = spreactral_density_lorentzian(fx,varargin)
	if (2 == length(varargin))
		f0 = varargin{1};
		p  = varargin{2};
	else
		f0 = varargin{1}(1);
		p  = varargin{1}(2);
	end

	S = 1./(1 + (p.*(abs(fx./f0)-1)).^2);
	% normalize
	if (isnumeric(fx))
	if (0)
		fmax = max(fx);	
		S = p./(f0*(atan(p*(fmax/f0 - 1)) + atan(p))).*S;
	else
		I = spectral_density_area(fx,S);
		S = S./I;
	end
	end
end
