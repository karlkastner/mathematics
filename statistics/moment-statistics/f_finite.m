% 2015-02-11 14:39:33.580789012 +0100
% Karl Kastner, Berlin
%
%% reduction factor of standard error for sampling from a finite distribution
%% without replacement
function f = f_finite(n_pop,n_sample)
%	f = sqrt((n_pop - n_sample)./(n_pop - 1));
	f = sqrt((n_pop - n_sample)./(n_pop.*n_sample));
end

