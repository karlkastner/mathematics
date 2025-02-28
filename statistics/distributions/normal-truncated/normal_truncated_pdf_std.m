% Tue 11 Feb 15:33:00 CET 2025
function sd = normal_truncated_pdf_std(mu,sd,lb,ub)
	mu = mu-lb;
	h = -mu/sd;

	r = normpdf(h)/(1-normcdf(h));
	% Z = normcdf(a)
	% s2  = 1 - a*normpdf(a)/Z - (normpdf(a)/Z).^2
	s2 = sd^2*(1 + r*(h - r)); 
	sd = sqrt(s2);
end

