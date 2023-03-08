% Wed 27 Apr 11:02:07 CEST 2022
function [S,R,r] = lowpass2d_pdf_discrete(L,n,La,p);
	% TODO this is only exact for odd
	fx = fourier_axis(L,n);
	x1  = fourier_axis(fx);
	x2  = x1;
	%x1 = ((1:n(1))-(n(1)+1)/2)'*dx;
	%x2 = ((1:n(2))-(n(2)+1)/2)*dx;
	r = hypot(cvec(x1),rvec(x2));
	% acf
	R = exp(-r./La);
	% spectal density
	%S = fft2(ifftshift(R));
	S = fft2(R);
	% should be real up to rounding error, just to make sure
	S = real(S);
	if (nargin()>3 && ~isempty(p))
		S = S.^p;
	end
	%S = S/sum(S(:)*dx^2);
	% normalize to 1 at 0
	S = S/S(1,1);
end


