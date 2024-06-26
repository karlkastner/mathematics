% Fri 17 May 08:48:42 CEST 2024
% Karl Kastner, Berlin
function p = rational_sigmoid_pdf(x,mu,a,k)
	x = abs(a.*(x-mu));
	p = 0.5*a./(1 + x.^k).^(1+1./k);
	% p = 0.5*((k-1)*a*x)/(k*(a*x + abs(a*x).^k).^(1 + 1/k);
end
