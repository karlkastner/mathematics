% Sun 11 Jul 12:19:45 CEST 2021
% c.f. umansky 2005, eq 13
% note that this is ident to:
% r(1)^2*D2x_c + r(1)*r(2)*(Dx_c*Dy_c+Dy_c*Dx_c) + r(2)^2*D2y_c;
function D = diffusion_matrix_2d_anisotropic(n,L,rho,a)
	if (length(n)==1)
		n = [n,n];
	end

	nn  = prod(n);
	h   = L./(n-1);
	r   = h(1)/h(2);

	buf = zeros(9*nn,3);

	c = cos(a);
	s = sqrt(1-c*c);

	A = c*c*rho(1) + s*s*rho(2); 
	B = s*s*rho(1) + c*c*rho(2);
	C = c*s*(rho(1)-rho(2));

	% diffusion stencil
	D = 1/h(1)^2*[-C*r/2,  B*r^2,       C*r/2
                     A,  -2*(A+B*r^2), A;
                     C*r/2,   B*r^2,      -C*r/2];
if (0)
	% this is a suggested symmetric higher order accurate alternative,
	% however, it is not fully clear to me why it works
	% the spectrum does also not reproduce the rotation better
	D = 1/h(1)^2*[0, B-C, C;
                   A-C, -2*(A+B-C), A-C;
	           C, B-C,0];
end	

	D = kernel2matrix(n,D);
end
