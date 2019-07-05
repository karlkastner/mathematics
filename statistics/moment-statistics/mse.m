% 2014-11-13 19:21:35.062315605 +0100
% Karl Kastner, Berlin
%
%% mean squared error of residual vector res
%% this is de-facto the std for an unbiased residual
function mse = mse(err,dim)
	% abs is necessary for imaginary values
	if (nargin<2)
		mse = nanmean(abs(err.*err));
	else
		mse = nanmean(abs(err.*err),dim);
	end
end

