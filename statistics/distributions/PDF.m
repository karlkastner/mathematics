% 2015-07-31 17:51:12.505066364 +0200
% Karl Kastner, Berlin
%
%% class for quasi-distributions from a set of sampling points
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
classdef PDF < handle
	properties
		n;
		CDF
		rank
		val
		s
		m
	end
	methods
	function obj = PDF(x)
		fdx = isfinite(x);
		obj.n = sum(fdx);
		[obj.val, obj.rank] = sort(x(fdx));
		obj.CDF = (1:obj.n)/(obj.n+1);
		% make values unique
		[obj.val udx] = unique(obj.val);
		obj.CDF = obj.CDF(udx);
		q = quantile(x,normcdf([-1 0 1]));
		obj.m = q(2);
		obj.s = 0.5*(q(3)-q(1));
	end
	function cdf = cdf(obj,x)
		cdf = interp1(obj.val,obj.CDF,x,'linear');
		% limit
		cdf = min(max(cdf,0),1);
	end
	function x = normalise(obj,x)
		 x = obj.s*norminv(obj.cdf(x)) + obj.m;
	end
	end %methods
end % class PDF

