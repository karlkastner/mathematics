% Fri  3 Dec 10:52:22 CET 2021
% Brownian noise
function [B,e,sS] = brownian_noise1d(L,n,e)
	if (length(n)==1)
		n(2) = 1;
	end
	if (nargin()<3)
		% uncorrelated Gaussian noise
		e = randn(n);
	end
	% padd inverse to ensure periodicity of integral
	e = [e;-e];
	dx = L/n(1);
	% angular frequency
	o = 2*pi*fourier_axis(L,n(1));
	% derivative operator
	D = 1i*o;
	% square root of spectral density of Brownian noise
	% (antiderivative operator)
	sS = 1./D;
	% set mean to zero
	sS(1) = 0;
	% Brownian noise
	% note, instead of fft(e), a complex random vector could be drawn
	B = 1/sqrt(dx)*ifft(sS.*fft(e));
	B = real(B);
	% extract first half
	B = B(1:n(1),:);
end

