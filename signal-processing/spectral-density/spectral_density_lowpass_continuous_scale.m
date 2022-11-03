% Tue 23 Aug 10:49:44 CEST 2022
function Sc = spectral_density_lowpass_continuous_scale(fc,p)
%	switch (p)
%	case {1}
%		Sc = 2./(pi*fc);
%	otherwise 
%		error('not yet implemented')	
		% I =  f*2F1(1/2,p,3/2,-f^2/fc^2);
		%   =  f*(1+f^2/fc^2)^(-1/2)*2F1(1/2,3/2-p,3/2,f^2/(f^2+fc^2)
		% at inf:
		% Iinf = fc*2F1(1/2,3/2-p,3/2)
		Sc = 1./(fc.*hypergeom([1/2,3/2-p],3/2,1));
%	end
end
