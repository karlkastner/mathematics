% Fri Feb 13 18:53:23 CET 2015
% Karl Kastner, Berlin
%% class for leave out 1 (delete 1) Jackknife estimates
%%
%% note 1 : the 1-delete jackknife does not yield consistend estimates for all functions,
%%        in particular it will perform poorly on robust estimation functions
%%        this is overcome by the d-delete jacknife, where d has to exceed the breakdown point
%%        of the estimating function, for example sqrt(n) for the median
%%        as this leads to unreasonably large number of repetitions, bootstrap
%%        is recommended for large sample cases (or blocking for sequential data)
%% note 2 : as a linearisation, jackknife underestimates the error variance in case of 
%%          dependence in the data
%% note 3 : studentisation and the leave out 1 jackknife are related
%% note 4 : the double 1 sample jacknife performs iferior to the d1 jacknife
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
classdef Jackknife < handle
	properties
		% function to jackknife, can return vectors
		% input must be organised in rows per sample
		% if actual function is multidimensional/indexed,
		% use a wrapper that has as first argument a mask
		func
		val
		bias
		serr
		C
		% R2 cannot be computed in general if the estimation function
		% has multiple multidimensional inputs acessed by masks/indices
		R2
		val0;
		val1;
		hat = [];
		hat0;
	end
	methods (Static)
		M = matrix1_STATIC(func,varargin);
		[param bias serr N C] = estimated_STATIC(param0,paramd,d,hat);
	end
	methods
		function obj = Jackknife(func)
			obj.func = func;
		end
		% estimates the parameters
		function obj = estimate(obj,varargin)
			% parameters with all data
			val0 = feval(obj.func,varargin{:});
			% parameters with leave out 1
			val1 = obj.matrix1(varargin{:});
			[val, bias, serr2 C] = obj.estimated_STATIC(val0,val1,1,obj.hat0);
			obj.val  = val;
			obj.bias = bias;
			obj.serr = sqrt(serr2);
			obj.val0 = val0;
			obj.val1 = val1;
			obj.C    = C;
%			obj.R2 = 1 - serr.*serr/
		end % estimate

		function [val, bias, res2, C, obj] = apply(obj,func,varargin)
			val0 = feval(func,obj.val0,varargin{:});
			for idx=1:size(obj.val1,3)
				val1(:,idx) = feval(func,obj.val1(:,:,idx),varargin{:});
			end
			if (nargout() > 3)
				[val, bias, res2 C] = obj.estimated_STATIC(val0, val1, 1, obj.hat0);
			else
				[val, bias, res2] = obj.estimated_STATIC(val0, val1, 1, obj.hat0);
			end
		end

		% jackknife matrix of pseudo variables
		function [M obj] = matrix1(obj,varargin)
			M = obj.matrix1_STATIC(obj.func,varargin{:});
		end % matrix1
	end % methods
end % classdef Jackknife

