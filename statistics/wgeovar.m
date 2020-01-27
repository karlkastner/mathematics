% Fri  3 Jan 15:45:06 +08 2020
% f = exp(sum wi log(xi)) = prod xi^wi
% s2 = (df/dxi)^2 e_i^2
% df/dxi = wi xi^(wi-1) f/(xi^w) = wi/xi f
% e_i = xi - f
%% variance of the weighted geometric mean
function s2 = wgeovar(w,x)
	w = w/sum(w);
	mu = wgeomean(w,x);
	% should the exponent here really be taken?
	%s2 = exp(sum(w.^2.*(log(x)-log(mu)).^2));
	%s2 = (sum((w-1)./(x.^(w-1)))).^2*mu.^2;
	%s2 = ((w-1).*(x.^(1-w))).^2*mu.^2.*(x-mu).^2;
	%s2 = sum((w./x).^2*mu.^2.*(x-mu).^2)
	s_div_mu2 = (w./x).*(x-mu);
	s2 = sum(s_div_mu2.^2).*mu.^2;
end
