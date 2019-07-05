% 2015-08-03 10:26:21.961541694 +0200
% Karl Kastner, Berlin
%
%% kurosis computed for quantiles
%%
%% Note : this is a measurement of shape-tailedness and yields the same value for the
%%        normal distribution as "kurtosis"
%%        However, this is a separate statistic and hence requires different
%%        methods for calculating P-values and hypothesis testing
function k = qkurtosis(x)
	q = quantile(x,normcdf([-1.5 -1 1 1.5]));
	if (isvector(q))
		q = cvec(q);
	end
	k = 2*(q(4,:)-q(1,:))./(q(3,:)-q(2,:));
end
