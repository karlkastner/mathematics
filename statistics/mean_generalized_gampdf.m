% Thu  2 Jul 17:31:40 +08 2020
% mean of the generalized gamma distribution
% E[X.^c], X := gamrnd(a,b)
% in particular for c = -1, the inverse gamma distribution
function mu = mean_generalized_gampdf(a,b,c)
	mu = b.^c.*gamma(a+c)./gamma(a);
end

