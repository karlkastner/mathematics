% Wed 10 Aug 16:15:53 CEST 2016
% Karl Kastner, Berlin
%
%% armijo stopping criterion for optimizations
%
% armijo 1966, eq 2
function sufficient = armijo_stopping_criterion(f0,g0,f,a)
	sufficient = false;
	% f(x - p*grad ) <= -1/2 a * g'*g
	if (f <= f0 - 0.5*a*(g0'*g0))
		sufficient = true;
	end
end

