% 2015-06-23 20:02:00.333016179 +0200
% Karl Kastner, Berlin
%
%% reduction of variance when sampling from a finite population
%% without replacement
function f = f_var_finite(n,m)
%	f = m.*(1./m - 1./n).^2 + (1./n).^2.*(n-m);
	f = (n-m)./(n.*m);
%	f = (n.*n - 3*n.*m + 2*m.*m)./(n.*n.*m);
end
