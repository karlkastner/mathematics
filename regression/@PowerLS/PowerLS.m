% Sat 15 Jul 16:18:41 CEST 2017
%% class for power law regression
classdef PowerLS < handle
	properties
		param
		rmse
		r2
		params
		nonlinear = true;
	end
	methods
		function obj = PowerLS()
		end
	end
end % PowerLS
