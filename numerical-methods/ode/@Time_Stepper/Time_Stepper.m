% Thu 20 Aug 10:50:22 +08 2020
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
classdef Time_Stepper < handle
	properties
		% start and end time of similartion (seconds) 
		Ti


		% numerical scheme for time stepping
		scheme = 'upwind';
		
		% number of terms for multi-step and runge-kutta
		order = 1;

		% factor for courant-friedrich-lewis stability condition
		cfl = 0.99;

		% interval at which the bed level is written to the output (seconds)
		% if 0, write always
		dto = 0;

		% stopping criteria for implicit schemes
		reltol = sqrt(eps);
	
		n_realloc = 100;
	end
	methods

	end % methods
end
