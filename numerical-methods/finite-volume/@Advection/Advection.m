% Tue 14 Nov 19:53:48 CET 2017
%
%% FVM treatment of the Advection equation
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

