% Wed 27 Apr 11:05:15 CEST 2022
function [S,R,r] = spectral_density_highpass2d_ideal(n,dx,L,p);
	[S,R,r] = spectral_density_lowpass2d_ideal(n,dx,L);
	q = 2;
	sS = S.^(1/q);
	sS = (1-sS);
	S  = sS.^q;
	if (nargin()>3)
		S = S.^p;
	end
end
