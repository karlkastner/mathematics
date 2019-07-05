% 2014-11-13 19:21:35.062315605 +0100
% Karl Kastner, Berlin
%
%% root mean square error computed from a residual vector
%% this is de-facto the std for an unbiased residual
function rmse = rmse(err,dim)
	if (nargin<2)
		rmse = sqrt(mean(err.*err));
	else
		rmse = sqrt(mean(err.*err,dim));
	end
end

