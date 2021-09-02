% 2021-06-27 01:09:43.796079281 +0200
function a = acf_bp(x,rho,dx,order)
	if (nargin()<3)
		dx = 1;
	end
	if (nargin()<4)
		order = 1;
	end
	switch (order)
	case {1}
	k   = x/dx;
	r   = 1/2*(1-rho^2)/rho;
	%r   = 2*(1-rho)/(1+rho);
	%a  := exp(-log(rho)*abs(k)).*(1-p*abs(k));
	a    = rho.^abs(k).*(1-r*abs(k));
	otherwise
		a = NaN(size(x));
	end
	% Discrete Communication Systems, Stevan Berber
	% only for sudden drop-off
%	a = sinc(p(1)*abs(x)).*cos(0.5*p(2)*abs(x));
end

