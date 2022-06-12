% Wed 27 Apr 11:02:07 CEST 2022
function [S,R,r] = spectral_density_lowpass2d_ideal(n,dx,L,p);
	% TODO this is only exact for odd
	x1 = ((1:n(1))-(n(1)+1)/2)'*dx;
	x2 = ((1:n(2))-(n(2)+1)/2)*dx;
	r = hypot(x1,x2);
	% acf
	R = exp(-r./L);
	% spectal density
	S = fft2(ifftshift(R));
	% should be real up to rounding error, just to make sure
	S = real(S);
	if (nargin()>3)
		S = S.^p;
	end
	%S = S/sum(S(:)*dx^2);
	% normalize to 1 at 0
	S = S/S(1,1);
end


