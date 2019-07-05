% 2014-11-13 19:21:35.062315605 +0100
% Karl Kastner, Berlin
%
%% root mean square error from vector of residuals
%% this is de-facto the std for an unbiased residual
function rmse = nanrmse(err,dim)
	if (nargin<2)
		rmse = sqrt(nanmean(err.*err));
	else
		rmse = sqrt(nanmean(err.*err,dim));
	end
end

