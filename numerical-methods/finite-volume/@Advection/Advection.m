% Tue 14 Nov 19:53:48 CET 2017
%
%% FVM treatment of the Advection equation
classdef Advection < handle
	properties (Constant)
		a = 1;
	end
	methods (Static)	
		
	function [dt] = dt_cfl(q,dx)
		dt = dx/a;
	end
	% roe average
	% TODO, make this a two-argument function ul and ur
	function [u] = roe_average(u)
		ul = [u(2:end);u(end)];
		u  = 0.5*(ul+u);
	end
	% jacobian matrix
	function [L R Ri] = fluxmateig(t,u)
		L  = Advection.a;
		R  = ones(size(u));
		Ri = ones(size(u));
	end
	% flux
	function [f] = flux(t,u)
		f = Advection.a*u;
	end
	function [A rhs] = bc_nonreflecting(t,q)
		A   = [0,1];
		rhs = 0;
	end
	end % methods (Static)
end % classdef Advection

