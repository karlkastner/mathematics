% 2016-10-23 20:54:10.409476622 +0200
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
%% transform modes (mu,sd) to parameters of the gamma distribution
%
function [a, b] = gampdf_moment2param(mu,sd)
	% mean     = a*b
	% variance = sd^2 = a*b^2
	b = sd.^2./mu;
	a = mu./b;
end

