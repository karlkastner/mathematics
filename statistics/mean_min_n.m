% Wed  1 Mar 17:18:17 CET 2023
function mean_min_n = mean_min_n(fun,Fun,n,varargin)
	mean_min_n = n*integral(@(x) x.*fun(x,varargin{:}).*(1-Fun(x,varargin{:})).^(n-1),0,inf);
end

