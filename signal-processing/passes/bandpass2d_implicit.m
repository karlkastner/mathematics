% Thu 24 Jun 14:39:10 CEST 2021
function [y] = bandpass2d_implicit(x,rho,a,order,direct)
	if (any(rho>=1))
		warning('rho must be smaller 1');
	end
	if (length(rho) < 2)
		rho = [rho,rho];
	end
	if (nargin()<5||isempty(direct))
		direct = false;
	end
	n  = size(x);
	tol = [];
	maxit = 2*max(n);
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,n-1,2,'circular');
	I = speye(prod(n));
	% y[i] = rho*(y[i-1] + y[i+1]) + (1-2*rho)*x[i]
	% (1 - 2 rho)y = rho*[1,-2,1]*y
	% ((1-4*rho)*I - rho*D2)*y = (1-4*rho)*x
	% y = (1-4*rho)*(((1-4*rho)*I - rho*D2))^-1*x
	% lowpass
	% A = (((1-2*(rhox+rhoy))*I - rhox*D2x - rhoy*D2y));

	% TODO is the normalization correct?
	r  = rho./(1 - 2*rho + rho.*rho)
	D2 = (r(1)*D2x + r(2)*D2y);

	if (nargin()>2&&~isempty(a))
		s = sin(a);
		c = cos(a);
		D2 = c*c*D2x + 2*s*c*Dxy + s*s*D2y;
		% TODO asymmetric scaling
		D2 = r(1)*D2;
	end

	A = (I - D2);
	y = flat(x);
	if (mod(order,1)==0)
	for idx=1:order
	if (direct)
		A2 = A*A;
		y = (A2 \ ((A-I)*(y)));
	else
		y =     pcg(A, y,tol,maxit);
		y = y - pcg(A, y,tol,maxit);
		if (0)
			% implementation in two steps is faster, better conditioned!
			A2 = A*A;
			y = pcg(A2,(A-I)*flat(y),tol,maxit,[],[],flat(x));
		end
	end
		%y = y/norm(y);
	end
	else
		A2 = A*A;
%		[V,E]=eigs(
	end
	y = reshape(y,n);
end

