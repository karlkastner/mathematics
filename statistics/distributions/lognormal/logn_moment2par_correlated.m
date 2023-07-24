% Wed 29 Mar 16:23:12 CEST 2023
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
function [lmu,lsd,lr] = logn_moment2par_correlated(mu,sd,r)
	[lmu,lsd]     = logn_moment2par(mu,sd);
	lr = log(1 + r.*(exp(lsd.^2) - 1))./lsd.^2;
end

