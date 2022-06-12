% Sun 11 Jul 20:28:23 CEST 2021
function [Dx,Dy,D2x,Dxy,D2y] = fourier_derivative_matrix_2d(n,L)
	s = -1;
	if (0 == mod(n(1),2))
		kx = 2i*pi*[0:n(1)/2-1,s*(n(1)/2:-1:1)]'/L(1);
	else
		kx = 2*pi*[0:(n(1)-1)/2,s*((n(1)-1)/2:-1:1)]/L(1);
	end
	if (0 == mod(n(2),2))
		ky = 2i*pi*[0:n(2)/2-1,s*(n(2)/2:-1:1)]'/L(2);
	else
		ky = 2*pi*[0:(n(2)-1)/2,s*((n(2)-1)/2:-1:1)]/L(2);
	end
	D1x = diag(sparse(kx));
	D1y = diag(sparse(ky));
	Ix  = speye(n(1));
	Iy  = speye(n(2));

	% y is first coordinate (columns) per matlab definition
	Dx = kron(Iy,D1x);
	Dy = kron(D1y,Ix);

	Dxy = kron(D1y,D1x);
	D2x = kron(Iy,D1x*D1x);
	D2y = kron(D1y*D1y,Ix);
end

