% Sun 12 Nov 18:33:42 CET 2017
%% class of flux limiters
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
classdef Flux_Limiter % not inherited from handle
	methods (Static)
		phi = beam_warming(theta);
		phi = fromm(theta);
		phi = lax_wendroff(theta);
		phi = minmod(theta);
		phi = monotized_central(theta);
		phi = superbee2(theta);
		phi = superbee(theta);
		phi = upwind(theta);
		phi = vanLeer(theta);
	end % static methods
	methods
		function obj = Flux_Limiter()
		end
	end
end % class Flux_Limiter

