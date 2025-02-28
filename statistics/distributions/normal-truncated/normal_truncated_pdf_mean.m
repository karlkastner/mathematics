% Tue 11 Feb 15:32:51 CET 2025
function mean_ = normal_truncated_pdf_mean(mu,sd,lb)
	error('broken')
	% shift lb to 0
	mu = mu-lb;
	h = -mu/sd;
	r = normpdf(h)/(1-normcdf(h));
	mean_ = mu + r*sd;
end

