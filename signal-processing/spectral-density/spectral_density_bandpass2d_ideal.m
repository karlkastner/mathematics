% Wed 27 Apr 11:05:35 CEST 2022
function [S,R,r] = spectral_density_bandpass2d_ideal(n,dx,Lf,p);
	[S,R,r] = spectral_density_lowpass2d_ideal(n,dx,Lf);
	q = 2;
%	sS = sqrt(S);
	sS = S.^(1/q);
	sS = sS.*(1-sS);
	S  = sS.^q;
	if (nargin()>3)
		S = S.^p;
	end
	% note that the integral for p = 1 is 2 pi L^2
	L_ = n.*dx;
	df = 1./L_;
	I = sum(sum(S))*df(1).^2
	S = 2*S./I;
end

