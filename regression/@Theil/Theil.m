% Fri Feb 20 11:11:51 CET 2015
% Karl Kastner, Berlin
%
%% Kendal-Theil-Sen robust regression
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
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

