% Fri 17 May 08:48:26 CEST 2024
% Karl Kastner, Berlin
function F = rational_sigmoid_cdf(x,mu,a,k)
	x = (a.*(x-mu));
	F = 0.5*(1+x./(1 + abs(x).^k).^(1./k));
end

