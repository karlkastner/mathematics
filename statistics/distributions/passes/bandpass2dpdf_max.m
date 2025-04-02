% Fri  3 Mar 13:08:40 CET 2023
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
%
%% Sc = max(S(fx,fy)) : maximum of the bandpass density in two dimensions in continuous space
%
function [Sc] = bandpass2dpdf_max(fc,p)
%	Sc = 1./(1/2*pi*fc^2*2^(4*p)/(2*p+1)*gamma(2*p+2)*gamma(2*p-1)/(gamma(4*p)));
	Sc = 2*(2*p+1)*gamma(4*p)./(pi*fc^2*2^(4*p)*gamma(2*p+2)*gamma(2*p-1));
end

