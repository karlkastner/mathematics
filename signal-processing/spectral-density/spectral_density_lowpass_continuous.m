% Sat 26 Jun 21:04:19 CEST 2021
function [S_lp, S_lp1] = spectral_density_lowpass_continuous(fx,fc,p,normalize)
	if (nargin()<3)
		p = 1;
	end
	if (nargin()<4)
		normalize = true;
	end
	S_lp1 = 1./(1 + fx.^2./fc.^2);
	S_lp  = S_lp1.^p;
	switch (normalize)
	case {0}
		% no normalization
		Sc = 1;
	otherwise
		Sc = spectral_density_lowpass_continuous_scale(fc, p);
	end
	S_lp = Sc.*S_lp;
end

