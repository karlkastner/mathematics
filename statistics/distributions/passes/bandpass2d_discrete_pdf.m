% Wed 27 Apr 11:05:35 CEST 2022
function [S,R,r] = bandpass2d_pdf_discrete(L,n,Lf,p,q);
	if (nargin()<5)
		q = 1;
	end
	
	[S,R,r] = lowpass2d_discrete_pdf(L,n,Lf);
	sS = S.^(1/q);
	% bandpass density
	sS = sS.*(1-sS);
	S  = sS.^q;
	% higher order
	if (nargin()>3 && ~isempty(p))
		S = S.^p;
	end
	% normalize
	% note that the integral for p = 1 is 2 pi L^2
	df = 1./L;
	I = sum(S,'all')*df(1)*df(2);
	S = 2*S./I;
end

