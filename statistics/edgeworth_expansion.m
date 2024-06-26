% Thu 17 Aug 12:24:20 CEST 2023
% Karl Kastner, Berlin
% note that the expansion is not guaranteed to be positivity preserving
function f = edgeworth_expansion(x,mu,sd,sk,ku,c5)
	if (nargin()<4)
		sk = [];
	end
	if (nargin()<5)
		ku = [];
	end
	n = 1;

	x_ = (x-mu)/sd;
	f = 1;
	m2 = sd^2;
	if (~isempty(sk))
		h3 = hermite(3,x_);
		m3 = sk*sd^(3);
		k3 = m3;
		l3 = k3/sd^3;

		f = f + 1/sqrt(n)*l3/6*h3;
	end
	if (~isempty(ku))
		m4 = ku*sd^4;
		k4 = m4 - 3*sd^4;
		l4 = k4/sd^4;
		h4 = hermite(4,x_);
		h6 = hermite(6,x_);
		f = f + 1/n*(l4/24*h4 + l3^2/72*h6);
	end
	if (~isempty(k5))
		h5 = hermite(5,x_);
		h7 = hermite(7,x_);
		h9 = hermite(9,x_);
		m5 = c5*sd^5;
		k5 = m5 - 10*m3*m2;
		l5 = k5/sd^5;
		f = f + 1/(n^(3/2))*(1/120*l5*h5 + 1/144*l3*l4*h7 + 1/1296*l3^3*h9);
	end
	f = f.*normpdf(x,mu,sd);
end

