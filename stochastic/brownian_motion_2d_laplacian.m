% 2022-06-13 10:18:31.647787389 +0200
%
% c.f. lodhia 2014
%
% H = 1/2
% s = H+d/2 <=> H = s-d/2
% Levy's Brownian motion: s = d/2 + 1/2 = 1 for 1d, 1.5 for 2d
% s_1d = 1
% s_2d = 1.5
%
% B = (-D)^s/2
%
% TODO : use analytical laplacian inverse with circular bc and truncate series

function [z,x1,x2] = brownian_motion_2d_laplacian(n,L)
	if (length(n)<3)
		n(3) = 1;
	end

	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,'circular');
	%[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,'circular');
	%[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,'dirichlet');
	%[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,L,2,'neumann');
	D2 = D2x+D2y;
	
	z = randn(n(1)*n(2),n(3));
	p = -1/2;
	p = -1;
	s = 1.5;
	p = s/2;
	if (1)
	%y =(-D2)^(p)*x; 
	z  =(-D2)^(p) \ z;
	z = z*3/2*sqrt(n(1)*n(2)); 
	%y =pcg(-D2,x); 
	else
		z = laplacian_power(L,n(1:2),-p,x);
	end
	x1 = (0:n(1)-1)'/n(1);
	x2 = (0:n(2)-1)'/n(2);
	%dx1 = x(2)-x(1);
	%dx2 = x2(2)-x2(1);

	z = reshape(z,n);

	% transform bridge into motion
	z = z + x1.*ones(1,n(2)).*randn(1,1,n(3));
	z = z + ones(n(1),1).*x2'.*randn(1,1,n(3));
	% mean
	z = z + randn(1,1,n(3));
end

