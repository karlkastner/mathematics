% 2023-03-30 09:50:37.717381573 +0200
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
function x_div_dx = int_1d_equal_exp(n,a_dx)
	% f = exp(-a*abs(x))
	% F = int_0^x f dx'/int_0^dx f dx' 
	%   = (1 - exp(-a*x))/(1 - exp(-a*dx))
	%   = (1 - exp(-a*x))/(1 - exp(-a*dx))
	%   = (1 - exp(-a_dx*x/dx))/(1 - exp(-a_dx))
	% qual inverval
	q = linspace(0,1,n+1);
	% q = F(x) -> x = F^-1(q)
	% Finv = log(1-q*(1-exp(-a_dx))/a_dx
	x_div_dx = log(1-q*(1-exp(-a_dx))/a_dx
end

