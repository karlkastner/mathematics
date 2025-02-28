% Thu  6 Jul 09:32:32 CEST 2023
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% initialize an empty object
function init(obj)
	obj.stat.area_msk = NaN;
	obj.stat.contrast = NaN;
	obj.stat.coverage = NaN;
	obj.stat.isisotropic = NaN;
	obj.stat.stati.intS_hp_sig = NaN;
	field_C = {'x','y','r'};
	for field = field_C
		obj.stat.L_eff.(field{1}) = NaN;
	end
	obj.stat.p_periodic = NaN;
	obj.stat.p_S_hp = NaN;
	field_C = {'hat','hp','bar','con'};
	for field = field_C
		obj.stat.fc.radial.(field{1}) = NaN;
		obj.stat.fc.x.(field{1}) = NaN;
		obj.stat.Sc.radial.(field{1}) = NaN;
		obj.stat.Sc.x.(field{1}) = NaN;
		obj.stat.Sc.y.(field{1}) = NaN;
	end
end

