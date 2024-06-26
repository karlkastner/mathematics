% Fri 23 Oct 15:51:08 +08 2020
% Karl Kastner, Berlin
function Phi = skewcdf(z,mu,sd,sk)
% T is owens t-functions (TFN)
%	Phi = normcdf(z) - 2*T(z,lambda);
	Phi = zeros(size(z),class(z));

	% asymptotic expansion : normcdf    ~ 1/2 + x/sqrt(2*pi)*exp(-x^2/2)
	% asymptotic expansion of : skewcdf = 2*sd_/sd*normcdf(a*x).*normpdf(x)	
	%                                   = 2*sd_/sd*(1/2+(a*x/(sqrt(2*pi))*exp(-(a*x)^2/2))*1/sqrt(2*pi)*exp(-x.^2/2)
	%                                   = 2*sd_/sd*(1/2*1/sqrt(2*pi)*exp(-x^2/2) + (a*x/((2*pi))*exp(-(a^2+1)*x^2/2))
	% note that exp(-x^2/2) ~ 1/(1+x^2/2) but this gets the tails wrong (overestimation)
	% the rhs (skew term) always decays faster
%	xl = fzero(@(x) skewpdf(x,mu,sd,sk)-sqrt(eps),mu-5*sd);
%	xr = fzero(@(x) skewpdf(x,mu,sd,sk)-sqrt(eps),mu+3*sd));
	xl = mu - 5*sd;	

	for idx=1:numel(z)
		xl_ = min(xl,z(idx));
		Phi(idx) = quad(@(x_) skewpdf(x_,mu,sd,sk),xl_,z(idx));
	end
end

