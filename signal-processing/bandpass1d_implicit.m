% Sat 26 Jun 10:35:59 CEST 2021
% TODO solve by cholesky factorisation / tridiagonal
function [y] = banspass1d_implicit(x,rho,order,invert)
	if (nargin()<4||isempty(invert))
		invert = false;
	end
if (0)
	% this is identical
	y = lowpass1d_implicit(x,rho,order,invert);
	y = y - lowpass1d_implicit(y,rho,order,invert);
else
	if (rho >= 1)
		warning('rho must be smaller 1');
	end
	if (nargin()<3)
		order = 1;
	end
	if (nargin()<4||isempty(invert))
		invert = false;
	end
	tol = [];
	n  = length(x);
	D2 = derivative_matrix_2_1d(n,n-1,[],'circular');         
	I  = speye(prod(n));
	% y[i] = rho*(y[i-1] + y[i+1]) + (1-2*rho)*x[i] 
	% (1 - 2 rho)y = rho*[1,-2,1]*y
	% ((1-2*rho)*I - rho*D2)*y = (1-2*rho)*x
	% y = (1-4*rho)*(((1-4*rho)*I - rho*D2))^-1*x
	% lowpass
	A = (I - rho/(1-2*rho+rho^2)*D2);
	maxit = max(n);

	y = x;
	for idx=1:order
	if (invert)
		A2 = A*A;
		y = (A2 \ ((A-I)*y));
	else
		if (1)
			y =     pcg(A,y,tol,maxit);
			y = y - pcg(A,y,tol,maxit);
		else
			% implementation in two steps is faster, better conditioned!
			A2 = A*A;
			y = pcg(A2,(A-I)*y,tol,maxit,[],[],y);
		end
	end
	end
end
end

