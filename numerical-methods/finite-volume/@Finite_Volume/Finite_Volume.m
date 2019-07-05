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

