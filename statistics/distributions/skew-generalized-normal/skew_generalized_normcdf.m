% Fri 23 Oct 15:51:08 +08 2020
% Karl Kastner, Berlin
function Phi = skew_generalized_normcdf(x,mu,sd,l1,l2)
% T is owens t-functions (TFN)
%	Phi = normcdf(z) - 2*T(z,lambda);
	Phi = zeros(size(x),class(x));

%	zl = fzero(@(z) skew_generalized_normpdf(z,mu,sd,l1,l2)-sqrt(eps),mu-5*sd);
%	xr = fzero(@(x) skewpdf(x,mu,sd,sk)-sqrt(e`ps),mu+3*sd));
	% for normal :     zl = 1/s*scale*exp(-(x-mu)/(2s))
	xl = mu-5*sd;

	for idx=1:numel(x)
		xl_ = min(xl,x(idx));
		Phi(idx) = quad(@(x_) skew_generalized_normpdf(x_,mu,sd,l1,l2),xl_,x(idx));
	end
end

