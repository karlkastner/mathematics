% Fri Feb 20 11:11:51 CET 2015
% Karl Kastner, Berlin
%
%% Kendal-Theil-Sen robust regression
%
classdef Theil < handle
	properties
		% intercep and slope
		param
		params
		% quantiles of X and Y and slope
		qx
		qy
		qi
		qs
		serr
		r2
		extended_statistics = true;
		repeated_medians = false;
		mad;
	end % properties
	methods (Static)
		%[slope, qs] = slope(X,Y,W,p);
		%[Y, slope]  = detrend(X,Y,W);
	end
	methods
		function obj = Theil()
		end
	end % methods
end % classdef Theil

