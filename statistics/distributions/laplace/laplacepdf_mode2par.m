% 2023-03-08 11:43:31.895412436 +0100
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
% function [mu,sd] = normpdf_mode2param(fc,Sc)
function [f0,s] = laplacepdf_mode2par(fc,Sc)
	if (issym(fc) || issym(Sc))
		pi_ = sym(pi);
	else
		pi_ = pi;
	end
	f0  = fc;
	s   = 1./(2*Sc);
	%mu = fc;
	%sd = 1./(Sc*sqrt(2*pi_));
end

