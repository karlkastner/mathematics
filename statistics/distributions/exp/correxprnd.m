% 2023-11-02 11:13:34.119189061 +0100
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% draw two vectors of correlated exponential random variables
% note: the solution is not unique(!)
% this solution is trivial, but not what one would like to have,
% var and rho are exact, but it is for every pair either 0, 1 not r
%
% note that x = p*x1 + sqrt(1-p^2)*x2 is not exponentially distributed 
%
function x = correxprnd(n,rho)
	if (rho<0)
		% correlations can only be positive in this sense
		error('impossible')
	end
	if (abs(rho)>1)
		error('impossible')
	end
	x = exprnd(1,n,2);
	id = randperm(n);
	id = id(1:round(rho*n));
	x(id,2) = x(id,1);
end

