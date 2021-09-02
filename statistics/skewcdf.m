% Fri 23 Oct 15:51:08 +08 2020
function Phi = skewcdf(z,mu,sd,sk)
% T is owens t-functions (TFN)
%	Phi = normcdf(z) - 2*T(z,lambda);
	Phi = zeros(size(z),class(z));
	

	xl = fzero(@(x) skewpdf(x,mu,sd,sk)-sqrt(eps),mu-5*sd);
%	xr = fzero(@(x) skewpdf(x,mu,sd,sk)-sqrt(eps),mu+3*sd));

	for idx=1:numel(z)
		xl_ = min(xl,z(idx));
		Phi(idx) = quad(@(x_) skewpdf(x_,mu,sd,sk),xl_,z(idx));
	end
end

