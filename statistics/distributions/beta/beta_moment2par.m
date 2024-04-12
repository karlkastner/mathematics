% 2016-10-23 21:12:46.030078012 +0200
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
%
%% transform central moments (mean and sd) to parameters of the beta function
function [a, b] = beta_moment_to_parameter(mu,sd)
	s2 = sd*sd;
	a = -mu*(mu*(1-mu)/s2-1)
	b = -(1 + (mu*mu - mu)/s2)*(mu - 1)

	a = abs(a);
	b = abs(b);

%	a = mu*(mu*(1-mu)/(s2)-1)
%	b = ((mu - 1)*(mu^2 + s2 - mu))/(s2)
%pause
%	b = (mu^3 - 2*mu^2 + mu*sd^2 + mu - sd^2)/sd^2
%	a =  b*mu/(1 - mu)
end

