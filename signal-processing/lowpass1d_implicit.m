% Sat 26 Jun 08:44:38 CEST 2021
% TODO solve by cholesky factorisation / tridiagonal
function [y] = lowpass1D_implicit(x,rho,order,invert)
	if (rho >= 1)
		warning('rho must be smaller 1');
	end
	if (nargin()<3||isempty(order))
		order = 1;
	end
	if (nargin()<4)
		invert = false;
	end
	n  = length(x);
	D2 = derivative_matrix_2_1d(n,n-1,[],'circular');         
	I  = speye(prod(n));
	% y[i] = rho*(y[i-1] + y[i+1]) + (1-2*rho)*x[i] 
	% (1 - 2 rho)y = rho*[1,-2,1]*y
	% ((1-4*rho)*I - rho*D2)*y = (1-4*rho)*x
	% y = (1-4*rho)*(((1-4*rho)*I - rho*D2))^-1*x
	% lowpass
	%A = (I - rho/(1-2*rho)*D2);
	A = (I - rho/(1-2*rho+rho^2)*D2);
	maxit = max(n);
	y = x;
	for idx=1:order
	if (invert)
		y = A \ y;
	else
		y = pcg(A,y,[],maxit);
	end
	end
end

