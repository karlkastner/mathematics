% Wed Apr 18 17:43:55 MSK 2012
% Karl Kästner, Berlin

% Derive FDM schemes as polynomials
% f_p = x_p * inv(X_sample)*f_sample ]
function D = derive_fdm_poly(x,x0)
	if (nargin() < 1)
		syms xk xl xc xr xs
		x = [xl; xc; xr]
	end
		
	if (nargin() < 2)
		syms x0
	end

	n = length(x);

	% vandermonde matrix at the grid points
	X = vander_1d(x,n-1);

	% n-th derivative of the vandermonde polynomial at x0
	syms B
	for idx=1:n
		% 1:n is obligatory, otherwise only first element is taken
		B(idx,1:n) = vanderd_1d(x0,idx-1,n-1);
	end

	% evaluate to get derivative coefficients at x0 (rows)
	%D = B / X % == X^-1*B
	D = B/X;
	%D = X\B
end % function derive_fdm_poly()

