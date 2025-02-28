% Tue 11 Feb 15:32:46 CET 2025
% this is truncated only at one end
function sk = normal_truncated_pdf_skewness(mu,sd,lb)
	error('broken')
	% shift lb to 0
	mu = mu+lb;
	h = -mu/sd;
	r  = normpdf(h)/(1-normcdf(h));
	sk = sd^3*r*((r-h).^2 + r*(r-h)-1);
end

