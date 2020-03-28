% Sun 19 Jan 22:07:40 +08 2020
% c.f. Wikipedia Mills Ratio
% E[X|X<Y] = mu - sigma phi(y)/Phi(Y)
function E = conditional_expectation_normal(y,mu,sigma)
	E = mu + sigma.*normpdf(y)./normcdf(y);
end

