% Wed 27 Apr 15:00:29 CEST 2022
function Sc = spectral_density_bandpass_2d_scale(fc,p,numeric)
	if (nargin()<3)
		numeric = false;
	end
	% normalize
	if (all(mod(p,1) == 0) && ~numeric)
		Sc = -3/(4*pi)*1./(16.^p.*fc.^2) ...
		     .* gamma(7/3-2*p)./(gamma(4/3-4*p).*gamma(2*p+1));
	else	
		% upper limit for integration
		frmax = cbrt(1e4*fc);
		fri = innerspace(0,frmax,1e5);
		dfr = fri(2)-fri(1);
		% radial integral
		Sc = 2*1./(2*pi*dfr*sum(fri.*spectral_density_bandpass_2d(fri,fc,p,false)));
	end
end

