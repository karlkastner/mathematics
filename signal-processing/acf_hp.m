% 2021-06-27 01:09:44.924080508 +0200
function a = acf_hp(x,rho,dx,order)
	if (nargin()<3)
		dx = 1;
	end
	if (nargin()<4)
		order = 1;
	end
	k = x/dx;
	% a := 1/2*(1-1/rho)*exp(-log(rho)*abs(x)/dx);
	switch (order)
	case {1}
		a  = 1/2*(1-1/rho)*rho.^abs(k);
		a(k == 0) = 1;
	case {2}
		%r   = 1/2*(1-rho^2)/rho;
		r  = 2*(1-rho)/(1+rho)
		a  = 1/2*(1-1/rho)*rho.^abs(k).*(1+r*abs(k));
		a(k == 0) = 1;
	otherwise
		a = NaN(size(x));
	end	
end

