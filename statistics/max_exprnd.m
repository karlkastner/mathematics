% 2020-04-28 13:20:27.003604475 +0800
% expected maximum of n iid exponential variables with mean mu
% E(max x) = sum 1/n ~ gamma + log(x) + 1/(2n)
function [m,mc] = max_exprnd(mu,n)
	m = mu*sum(1./(1:n));
	if (nargout()>1)
		mc = mu*cumsum(1./(1:n));
	end
end

