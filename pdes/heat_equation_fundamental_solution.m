% Wed 14 Jun 14:56:42 CEST 2023
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
%%
%% function y = heat_equation_fundamental_solution(t,x,d,t0,x0)
function y = heat_equation_fundamental_solution(t,x,d,t0,x0)
	y = normpdf(x,x0,sqrt(2*d*(t+t0)));
end

