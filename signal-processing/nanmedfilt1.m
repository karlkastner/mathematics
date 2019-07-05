% 2018-05-28 16:26:06.194124055 +0200 nanmedfilt1.m
%% medfilt1, skipping nans
function x = nanmedfilt1(x,varargin)
	fdx = 
	x = isfinite(x);
	x = medfilt1(double(x))'
end

