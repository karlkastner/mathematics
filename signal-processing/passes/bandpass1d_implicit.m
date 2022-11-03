% Sat 26 Jun 10:35:59 CEST 2021
% Karl Kastner, Berlin
%
% filter the input vector with a spatial (two-sided) bandpass
% in real space (finite difference approximation)
%
% function [y] = banspass1d_implicit(x,rho,order,invert)
% TODO solve by cholesky factorisation / tridiagonal
function [y] = banspass1d_implicit(x,rho,order,invert)
	if (nargin()<4||isempty(invert))
		invert = false;
	end
if (0)
	% this is identical
	y =     lowpass1d_implicit(x,rho,order,invert);
	y = y - lowpass1d_implicit(y,rho,order,invert);
else
	if (any(rho >= 1))
		warning('rho must be smaller 1');
	end
	if (nargin()<3||isempty(order))
		order = 1;
	end
	if (1 == length(rho))
		rho = [rho,rho];
	end
	if (isvector(x))
		x = cvec(x);
	end
	n  = length(x);
%	D2 = derivative_matrix_2_1d(n,n-1,[],'circular');         
	Dl = derivative_matrix_1_1d(n,n-1,-1,'circular');         
	Dr = derivative_matrix_1_1d(n,n-1,+1,'circular');         
	I  = speye(prod(n));
	% y[i] = rho*(y[i-1] + y[i+1]) + (1-2*rho)*x[i] 
	% (1 - 2 rho)y = rho*[1,-2,1]*y
	% ((1-2*rho)*I - rho*D2)*y = (1-2*rho)*x
	% y = (1-4*rho)*(((1-4*rho)*I - rho*D2))^-1*x
	% lowpass
	%r = rho/(1-2*rho+rho^2);
	r = rho./(1-rho(1)-rho(2)+rho(1)*rho(2));
	% note: D2 = Dr - Dl = D2
	A = (I - r(1)*Dl + r(2)*Dr);
	tol = [];
	maxit = max(n);
	y = x;
	if (abs(rem(order,1))>sqrt(eps))
		%A2 = A*A;
		%[v,e] = eig(A2);
		%y = (A2^order \ ((A-I)^order*y));
		warning('not yet implemeneted')
	end
	for idx=1:order
	if (invert)
		% TODO solve by cholesky factorisation
		A2 = A*A;
		y = (A2 \ ((A-I)*y));
	else
		if (1)
			% implementation in two steps is faster, better conditioned!
			y =     pcg(A,y,tol,maxit);
			y = y - pcg(A,y,tol,maxit);
		else
			% TODO precondition
			A2 = A*A;
			y = pcg(A2,(A-I)*y,tol,maxit,[],[],y);
		end
	end % if invert
	end % for idx
end % if 0
end % bandpass1d_implicit

