% 2012 Apr 27 15:49 (MSK)
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
% demonstration, that averaging densities with the same distribution
% but different regularity results in a density that is more pointed
% and has heavier tales than the underlying distribution

% coordinates and weights for numerical Gauss quadrature
% 3rd order accurate
function [w, b, flag] = int_1d_gauss_2()

	% quadrature points in baricentric coordinates	
	%a = 0.5*sqrt(1./3) + 0.5
	%b = [   a (1-a);
        %    (1-a)   a];
	% w = [1];
	% b = [0.577350269189626];
	% weights
 	w = [ 0.500000000000000; 0.500000000000000];
	b = [ 0.211324865405187 0.788675134594813;
              0.788675134594813 0.211324865405187];

	% no diagonal mass-matrix
	flag = 0;
end % int_1d_gauss_2()

