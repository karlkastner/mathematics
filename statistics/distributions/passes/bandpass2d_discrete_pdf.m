% Wed 27 Apr 11:05:35 CEST 2022
function [S,R,r] = bandpass2d_pdf_discrete(n,dx,Lf,p,q);
	if (nargin()<5)
		q = 1;
	end
	[S,R,r] = lowpass2d_discrete_pdf(n,dx,Lf);
	sS = S.^(1/q);
	sS = sS.*(1-sS);
	S  = sS.^q;
	if (nargin()>3)
		S = S.^p;
	end
	% note that the integral for p = 1 is 2 pi L^2
	L_ = n.*dx;
	df = 1./L_;
	I = sum(S,'all')*df(1).^2
	S = 2*S./I;
end

