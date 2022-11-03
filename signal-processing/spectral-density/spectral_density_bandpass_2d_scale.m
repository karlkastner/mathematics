% Mon 22 Aug 10:53:04 CEST 2022
% scale (normalization) factor of the spectral density of the 2D bandpass filter
% the scale factor is equal to the maximum of the density

% Wed 27 Apr 15:00:29 CEST 2022
function Sc = spectral_density_bandpass_2d_scale(fc,p,pp,numeric)
	if (nargin()<4)
		numeric = false;
	end
	if (nargin()<3)
		pp = [];
	end
	if (isempty(pp) && ~numeric)
		kc = 2*pi*fc;
		IS  = kc.^2./(2*pi) .* 2.^(2*p-1)./(p-1)./binom(2*p-1,p-1); 
		Sc = 1./IS;
	else
		error('not yet implemented')
	end


%	if (all(mod(p,1) == 0) && ~numeric)
%		Sc = -3/(4*pi)*1./(16.^p.*fc.^2) ...
%		     .* gamma(7/3-2*p)./(gamma(4/3-4*p).*gamma(2*p+1));
%	else	
%		% upper limit for integration
%		frmax = cbrt(1e4*fc);
%		fri = innerspace(0,frmax,1e5);
%		dfr = fri(2)-fri(1);
%		% radial integral
%		Sc = 2*1./(2*pi*dfr*sum(fri.*spectral_density_bandpass_2d(fri,fc,p,false)));
%	end
end

%function Sc = spectral_density_bandpass_2d_scale(fc,p)
%	kc = 2*pi*fc;
%	if (0)
%		I = kc^2*2^p./(p+2).*gamma(p/2 + 2).*gamma(p/2 - 1) ./ gamma(p);
%		Sc = 1./I;
%	else
%		Sc = (p+2)./2.^p./kc.^2.*gamma(p)./gamma(p/2-1)./gamma(p/2+2);
%	end
%end
%
