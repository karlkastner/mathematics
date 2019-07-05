% Sun 30 Jul 15:33:55 CEST 2017
%
%% expected maximum of n normal variables
%% c.f. Wolperts
%% this is the median, not the mean of the maximum!
%% see median of gumbel
function m = maxnnormal(mu,sigma,n)
	m = mu - sigma*norminv(1/n) + log(log(2))/norminv(1/n);
end

