% Sat Feb 19 19:59:39 MSK 2011
% Thu  2 Nov 21:01:36 CET 2017 (change to object)
% Karl Kästner, KTH Stockholm
%
%% finite volume method for partial differential equations 1+1 dimensions
%% (time and space)
%%
%
% conservation law -> balance law -> sink term, 19.5 splitting for 2d, 17.14
% TODO: for stationary flow solver:
%	build highres as matrix then solve
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
classdef Finite_Volume < handle
	properties
		step
		pde
%		flux
%		fluxmat

		% number of unknowns (1 for advection, 2 for swe, 3 for gas dynamics)
		m = 2

		dx
		x

		icfun	    = @(x) zeros(2*length(x),1);
		bcfun	    = {};
		%sourcefun   = @(t,x,q) zeros(2*length(x),1);
		cflscale    = 0.99;

		dt_progress = 10;

		odestepper  = @step_trapezoidal;

		startup     = true;

		dt_out      = 0;

		icabstol    = 1e-7;
		icreltol    = 1.e-6;
		%zbfun;
	end % properties

	methods (Static)
	
	end % static methods

	methods
		function obj = init(obj,X,n)
			L = diff(X);
			obj.dx = L/(n-1);
			obj.x = X(1) + obj.dx*(0:n-1)';
 			obj.step = @obj.step_unsplit;
			obj.pde.init(obj.x);
			% split_strang;
		end
	end % methods
end % class Finite_Volume

