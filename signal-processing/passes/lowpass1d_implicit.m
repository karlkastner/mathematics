% Sat 26 Jun 08:44:38 CEST 2021
% function [y] = lowpass1D_implicit(x,rho,order,invert)
% note : this is a script for demonstration, for efficient execution,
% use fourier-implementation
% function [y] = lowpass1D_implicit(x,rho,order,invert)
function [y] = lowpass1D_implicit(x,rho,order,invert)
	if (nargin()<4 || isempty(invert))
		invert = false;
	end

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
	% lowpass:
	% y[i] = rho*(y[i-1] + y[i+1]) + (1-2*rho)*x[i] 
	% (1 - 2 rho)y = rho*[1,-2,1]*y + (1-2*rho)*x[i]
	% ((1 - 2*rho)*I - rho*D2)*y = (1-2*rho)*x
	% y = (1-2*rho)*(((1-2*rho)*I - rho*D2))^-1*x
	% y = ((I - rho/(1-2*rho)*D2))^-1*x

% must hold for impulse response: x = [0...,1,...0], y = rho^[|k|]
% (I - r*D2)*y = x
% A*y = x
%(1 - r*[1,-2,1])*1/s*[rho,1,rho] = 1
%(1 - 2*r*(rho-1)) == s
%s = sum_-inf^inf r^|k| = 1/(1-rho) + rho/(1-rho)
%(1 - 2*r*(rho-1)) == 1/(1-rho) + rho/(1-rho)
% - 2*r*(rho-1) == -1 + 1/(1-rho) + rho/(1-rho)
% - 2*r == -1/(rho-1) - 1/(1-rho)^2 - rho/(1-rho)^2
% - 2*r ==  1/(1-rho) - 1/(1-rho)^2 - rho/(1-rho)^2
% - 2*r ==  (1-rho)/(1-rho)^2 - 1/(1-rho)^2 - rho/(1-rho)^2
% - 2*r ==  ((1-rho) - 1 - rho)(1-rho)^2
% - 2*r ==  (-2 rho)(1-rho)^2
%     r ==  rho/(1-rho)^2
%(1 - 2*r*(rho-1)) == s
%r = 1/2*(1-s)/(rho-1)
%
% 2D:
%(1 - r*[1,..,1,-4,1,...,1])*1/s*[rho,...,rho,1,rho,...,rho] = 1
%(1 - 4*r*(rho-1)) == s
% r = 1/4*(1-s)/(rho-1)

% y = rl yl + rr yr + e
% (1 - 1/2 rl + 1/2 rr) y = 1/2 rl (yl - yr) + 1/2 rr(yr - y) + e
% (1 - 1/2 rl + 1/2 rr) y = 1/2 rl (yl - yr) + 1/2 rr(yr - y) + e
% (1 - 1/2 rl + 1/2 rr) y = 1/2 rl (yl - yr) + 1/2 rr(yr - y) + e
% (I - (1/2 rl Dl + 1/2 rr Dr)/(1 - 1/2 rl + 1/2 rr)) y = e

	% note : this preserves the sum of x (l1-norm when strictly positive),
	%        not the rms!
	r = 0.5*rho./(1 - 0.5*rho(1) - 0.5*rho(2));
	% note: D2 = Dr - Dl = D2
	A = (I - r(1)*Dl + r(2)*Dr);
	tol = [];
	maxit = max(n);
	y = x;
	% note that this only applies to integer order
	for idx=1:order
	if (invert)
		% TODO solve by cholesky factorisation
		y = A \ y;
	else
		% TODO precondition
		y = pcg(A,y,tol,maxit,[],[],y);
	end
	end
end

