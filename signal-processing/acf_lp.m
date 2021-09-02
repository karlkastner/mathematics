% 2021-06-27 01:09:46.000081606 +0200
function a = acf_lp(x,rho,dx,order)
	if (nargin()<3)
		dx = 1;
	end
	if (nargin()<4)
		order = 1;
	end
	k = x/dx;
	%a = exp(log(rho)*abs(x)/dx)
	switch (order)
	case {1}
		a = rho.^abs(k);
	case {2}
		%K = max(k);
		%r   = 1/2*(1-rho^2)/rho;
		r   = 2*(1-rho)/(1+rho)
		a = rho.^abs(k).*(1+r*abs(k))
	otherwise
		a = NaN(size(x));
	end
end

