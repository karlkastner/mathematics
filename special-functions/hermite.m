% Do 10. Mär 14:29:35 CET 2016
% Karl Kastner, Berlin
%
%% probabilistic's hermite polynomial by recurrence relation
%%
%% input :
%% n : order
%% x : value
%%
%% output:
%% f : H_n(x)
%% df : d/dx H_n(x)
%%
function [f, df] = hermite(n,x)
	[f, df] = hermite_(n);
	if (nargin() > 1)
		f  = polyval_(f,x);
		df = polyval_(df,x);
	end
end

function f = polyval_(p,x)
	f = 0;
	n = length(p);
	for idx=1:n
		f = f + x.^(n-idx)*p(idx);
	end
end

function [f, df] = hermite_(n)
	if (n <= 0)
		f  = 1;
		df = 0;
	else
		[f_ df_] = hermite_(n-1);
		f  = [f_, 0] - [0, df_];
		df = [0, f(1:end-1).*(n:-1:1)];
	end
end


