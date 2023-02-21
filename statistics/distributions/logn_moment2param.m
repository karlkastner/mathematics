% Mi 2. Mär 14:47:56 CET 2016
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
% demonstration, that averaging densities with the same distribution
% but different regularity results in a density that is more pointed
% and has heavier tales than the underlying distribution
%
%% transform the mode (mu,sd) to parameters of the log normal distribution
%
function [lmu,lsd] = logn_moment2param(mu,sd)
	lsd = sqrt(log(1 + sd.^2./mu.^2));
	lmu = log(mu)-0.5*lsd.^2;
end
