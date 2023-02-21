% Fri  7 Jan 12:44:47 CET 2022
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
%% transform mode to parameters of the brownian phase spectral density
%
% function [f0, s] = spectral_density_brownian_phase_mode2par(fc,Sc)
function [f0, s] = spectral_density_brownian_phase_mode2par(fc,Sc)

	%r = root(z^3 - z^2 + z*(4*Sc^2*fc^2*pi^2 - 2) - 4*Sc^2*fc^2*pi^2, z, 1)
	rp = [1, -1, (4*Sc^2*fc^2*pi^2 - 2), -4*Sc^2*fc^2*pi^2];
	r = roots3(rp);
	% choose real root
	r = r(1);
	p2 = (  pi^2*((4*Sc^2*fc^2*pi^2 - 1)*r^2 ...
	      - 12*Sc^2*fc^2*pi^2 + 16*Sc^4*fc^4*pi^4 ...
              + 2*r)/(16*Sc^2*fc^2) );
	s = sqrt(pi/sqrt(p2));
	f0 = fc./sqrt(2*(pi.^2*s.^4 + 1).^(1/2) - s.^4*pi^2 - 1);
end % spectral_density_brownian_phase_mode2par

