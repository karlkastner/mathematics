% Fri  3 Dec 10:52:22 CET 2021
%
% generate a Brownian surface (Gaussian field) in two dimensions
%
% BM = Delta^(-1.5/2) e, where e is iid Gaussian
%    = F^-1 S F e, where S = |omega|^-(1.5/2)
%
% since it is expensive to compute the eigenpairs of the laplacian Delta,
% the fourier transform is employed
%
% the fourier transform yields a brownian bridge, due to the implied periodicity
%
function [B,e,sS] = brownian_noise_2d_fft(n,L,bridge,e)
	if (length(n)<3)
		n(3) = 1;
	end
	if (nargin()<2)
		L = [1,1];
	end
	if (nargin()<3)
		bridge = false;
	end
	if (nargin()<4)
		% uncorrelated Gaussian noise
		% we do not need to take the fourier transform of e,
		% as the fourier transform of gaussian noise just yields
		% another gaussian random vector
		e = randn(n(1),n(2),n(3));
	end
	% angular frequency
	% is this not wrong, bc: (Fx^-1 D2x Fx + Fy^-1 D1y Fy)^-1 is not F2^-1 (D2x + D2y) F2^-1 ?
	ox = 2i*pi*fourier_axis(L(1),n(1));
	oy = 2i*pi*fourier_axis(L(2),n(2));

	% radial angular frequency squared
	or2  = (cvec(ox).^2+rvec(oy).^2);	
	sS  = or2.^(-1.5/2);
	% the mean is an arbitrary integration constant, set to zero
	sS(or2 == 0) = 0;

	% fourier transform of brownian motion
	B = (sS.*e);
	% brownian motion
	B = real(ifft2(B));

	% TODO what is the exact scale factor to obtain the desired std?
	B = 5/2*(n(1)*n(2)).^(1)*real(B)*(n(1)*n(2)/((n(1)+1)*(n(2)+1)))^-(0.75/2);

	% transform bridge into motion
	if (~bridge)
		x = (0:n(1)-1)'/n(1);
		y = (0:n(2)-1)/n(2);
		B = B + sqrt(L(1))*x.*ones(1,n(2)).*randn(1,1,n(3));
		B = B + sqrt(L(2))*ones(n(1),1).*y.*randn(1,1,n(3));
	end

end % brownian_noise_2d_fft

