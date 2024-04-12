% Wed  9 Nov 12:35:34 CET 2022
% Karl Kästner, Berlin
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%% function sy = spectral_density_brownian_phase_across_mode2par(Scy)
%% by definition, the maximum occurrs at the oricin fc_y = 0
%
function sy = phase_drift_parallel_pdf_mode2par(Scy)
	sy = sqrt(4./Scy);
end

