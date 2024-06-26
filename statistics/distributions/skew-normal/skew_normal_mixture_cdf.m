% Thu 16 May 00:25:08 CEST 2024
% Karl Kastner, Berlin
function Phi = skew_normal_mixture_cdf(z,p,mu,sd,sk)
% T is owens t-functions (TFN)
%	Phi = normcdf(z) - 2*T(z,lambda);
	Phi = zeros(size(z),class(z));

	xl = min(mu-5*sd);

%	xl = fzero(@(x) skewpdf(x,min(mu),max(sd),sk)-sqrt(eps),mu-5*sd);
%	xr = fzero(@(x) skewpdf(x,mu,sd,sk)-sqrt(eps),mu+3*sd));

	p = p/sum(p);
	for idx=1:numel(z)
		xl_ = min(xl,z(idx));
		Phi(idx) = quad(@(x_) skew_normal_mixture_pdf(x_,p,mu,sd,sk),xl_,z(idx));
	end
end

