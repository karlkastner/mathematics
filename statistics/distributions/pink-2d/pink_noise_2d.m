% Mon 22 Jan 09:11:43 CET 2024
% Karl Kastner, Berlin
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
%
% noise where spectral energy decays by one order per octave
% note that this is different to brownian noise (see brownian noise generator)
%
%%% function [e,T,fx,fy,frr] = pink_noise_2d(n,L);
function [e,T,fx,fy,frr] = pink_noise_2d(n,L);
	[e,T,fx,fy,frr] = coloured_noise_2d(n,L,2.0);
end

