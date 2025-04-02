function [a,b,p] = generalized_gamma_cmoment2par(mu,s2,sk)
	% initial condition
	[a,b] = gamma_moment2par(mu,s2);
	par0 = [a,b,1];
	par = lsqnonlin(@(par) [generalized_gamma_mean(p(1),p(2),p(3)) - mu,
			  generalized_gamma_var(p(1),p(2),p(3)) - s2,
			  generalized_gamma_skewnewss(p(1),p(2),p(3)) - sk], ...
			  par0 ...
		 );
	a = par(1);
	b = par(2);
	p = par(3);
end

