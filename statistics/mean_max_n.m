% Wed  1 Mar 17:19:32 CET 2023
function mean_max_n = mean_max_n(fun,Fun,n,varargin)
	mean_max_n = n*integral(@(x) x.*fun(x,varargin{:}).*Fun(x,varargin{:}).^(n-1),0,inf);
end

