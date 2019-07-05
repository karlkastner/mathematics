% 2018-02-02 13:39:03.497357890 +0100
%% effective sample size correction for autocorrelated series
function n = ar1_effective_sample_size(rho,twosided)
	if (~twosided)
		n = (1+rho)/(1-rho);
	else
		n = (1 + rho)^3/((1-rho)*(1+rho^2));
	end
end

