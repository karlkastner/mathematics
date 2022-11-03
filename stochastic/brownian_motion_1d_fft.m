% Fri  3 Dec 10:52:22 CET 2021
% Brownian noise
function [B,e,sS] = brownian_noise_1d_fft(L,n,e)
	if (length(n)==1)
		n(2) = 1;
	end
	if (nargin()<3)
		% uncorrelated Gaussian noise
		e = randn(2*n(1),n(2));
	end
	% angular frequency
	o = 2i*pi*fourier_axis(2*L,2*n(1));
	% second derivative operator
	D2 = (o).^2;
	% square root of spectral density of Brownian noise
	% (antiderivative operator)
	sS = 1./sqrt(D2);
	% set mean to zero
	sS(1) = 0;
	% Brownian noise
	% note: 2*sqrt(n)/dx = 2*n/L
	B = 2*n(1)/sqrt(L)*ifft(sS.*e);
	B = real(B);
	% extract first half
	B = B(1:n(1),:);
if (1)
	% transform the brownian bridge into a brownian motion
	r = randn(1,n(2));
	x = sqrt(L)*(0:n(1)-1)'/(n(1)-1/2);
	B = B+x*r;
end
end

