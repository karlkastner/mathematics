% Wed Apr 18 17:43:55 MSK 2012
% Karl Kästner, Berlin

% 5-point FDM schemes by polynomials
% f_p = x_p * inv(X_sample)*f_sample ]
function D = derive_poly(x,x0)
	if (nargin > 0 && ~isempty(x))
		if (length(x) == 3)
			xl = x(1);
			xc = x(2);
			xr = x(3);
		else
			xk=x(1);
			xl=x(2);
			xc=x(3);
			xr=x(4);
			xs=x(5);
		end
		
	else
		syms xk xl xc xr xs x0
	end
	
	if (length(x) == 3)
	X = [1 1 1 
	     xl xc xr 
	     xl.^2 xc.^2 xr.^2 ].';
        B = [x0^0  x0^1   x0^2;
                0  x0^0 2*x0^1;
                0     0 2*x0^0];
	else
	X = [1 1 1 1 1
	 xk xl xc xr xs
	 xk.^2 xl.^2 xc.^2 xr.^2 xs.^2
	 xk.^3 xl.^3 xc.^3 xr.^3 xs.^3
	 xk.^4 xl.^4 xc.^4 xr.^4 xs.^4 ].';
	
	% evaluate n-th derivatives of the polynomials
	B = [x0^0  x0^1   x0^2   x0^3    x0^4; % yields identity
	        0  x0^0 2*x0^1 3*x0^2  4*x0^3; % yields first derivative
	        0     0 2*x0^0 6*x0^1 12*x0^2;
	        0     0      0 6*x0^0 24*x0^1;
	        0     0      0      0 24*x0^0 ];
	end

	% get the derivatives (rows)
	D = B / X; % == X^-1*B
end % function poly

