% Sun 11 Jul 20:28:23 CEST 2021
% TODO dx
function [Dp] = fourier_derivative_matrix_1d(n,L,p)
	if (nargin()<3)
		p = 1;
	end
	% TODO odd
	s = -1;
	kx = 2i*pi*[0:n(1)/2-1,s*(n(1)/2:-1:1)]'/(L(1));

	Dp = diag(sparse(kx.^p));
end

