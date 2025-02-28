% Wed 29 Mar 16:59:42 CEST 2023
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
function sd = geometric_ou_2d_grid_cell_averaged_std(lmu,lsd,theta,dx,dy,varargin)
	%lrfun = @(x,y) exp(-hypot(x,y)/theta);
	%% mean of the values val
	%mu = exp(lmu+0.5*lsd^2);
	%% covariance fucntion of the values
	%cfun = @(x,y) logn_corr(lrfun(x,y),lmu,lmu,lsd,lsd);
	s2 = geometric_ou_2d_grid_cell_averaged_cov(lmu,lsd,theta,0,0,dx,dy,varargin{:});
	sd = sqrt(s2);
end

