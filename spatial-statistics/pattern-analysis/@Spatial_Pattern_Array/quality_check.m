% 2022-09-26 14:32:10.449630621 +0200
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
function [qc_passed, qc_flags] = quality_check(obj)

	switch (obj.type)
	case {'isotropic'}
		wavelength = obj.wavelength_r;
		regularity = obj.regularity_r;
		isisotropic = 1;
	case {'anisotropic'}
		wavelength = obj.wavelength_x;
		regularity = obj.regularity_x;
		isisotropic = 0;
	end % switch type

	qc_flags = [    cvec(isfinite(cvec(obj.p_periodic))) ...
		      , cvec(obj.isisotropic == isisotropic) ...
		      , cvec(regularity >= obj.opt.regularity_min) ...
		      , cvec(isfinite(wavelength)) ...
		      , (cvec(wavelength)<obj.opt.wavelength_max) ...
		      , (cvec(wavelength)>obj.opt.wavelength_min) ...
		    ];

	qc_passed = all(qc_flags,2);
end % quality_check

