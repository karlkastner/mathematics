% Sun 11 Jul 20:28:23 CEST 2021
%
% function [Dp] = fourier_derivative_matrix_1d(n,L,p)
function [Dp] = fourier_derivative_matrix_1d(n,L,p)
	if (nargin()<3)
		p = 1;
	end
	s = -1;
	if (0 == mod(n,2))
		kx = 2i*pi*[0:n/2-1,s*(n/2:-1:1)]'/L(1);
	else
		kx = 2*pi*[0:(n-1)/2,s*((n-1)/2:-1:1)]/L(1);
	end

	Dp = diag(sparse(kx.^p));
end

