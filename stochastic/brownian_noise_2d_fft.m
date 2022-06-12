% Fri  3 Dec 10:52:22 CET 2021
% 2022-01-26 13:21:23.510203267 +0100
% Brownian noise
%
% dy/dx = d/dx F^-1 sum y_i exp(i omega x) =  sum i omega x 
% <=> int f dx = int F^-1 y_i exp(i omega x) = F^-1 sum 1/(i omega_i) y_i
%
%function [B,e] = brownian_noise_2d_fft(L,n,e)
function [B,e] = brownian_noise_2d_fft(L,n,e)
	if (length(L)==1)
		L(2) = L(1);
	end
	if (length(n)==1)
		n(2) = n(1);
	end
	if (nargin()<3)
		% uncorrelated Gaussian noise
		e = randn(n);
	end
	% padd iverse in both dimensions to ensure periodicity of integral
	% note, this works also for 3D matrices
	e = [e, -e];
	e = [e; -e];
	dx = L(1:2)./n(1:2);
	% angular frequency
	ox = 2*pi*fourier_axis(2*L(1),2*n(1));
	oy = 2*pi*fourier_axis(2*L(2),2*n(2));
	% radial frequency
	fr = hypot(cvec(ox),rvec(oy));
	% square root of spectral density of 2D brownian noise
	% (antiderivative operator)
	sS = 1./fr;
	% set mean to zero
	sS(1,1) = 0;
	% Brownian noise
	% note, instead of fft2(e), a complex random vector could be drawn
	B = 1./sqrt(dx(1)*dx(2))*ifft2(sS.*fft2(e));
	B = real(B);
	% extract the first quadrant
	B = B(1:n(1),1:n(2),:);
end

