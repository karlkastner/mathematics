% Wed 27 Apr 11:02:07 CEST 2022
function [S,R,x,y] = lowpass2d_pdf_discrete(L,n,La,p);
	% autocorrelation
	[R,x,y] = lowpass2d_discrete_acf(L,n,La);

	% spectal density
	S = fft2(R);

	% should be real up to rounding error, just to make sure
	S = real(S);

	if (nargin()>3 && ~isempty(p))
		S = S.^p;
	end
	% normalize to 1 at 0
	S = S/S(1,1);
end

