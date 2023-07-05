% Wed 29 Mar 16:58:33 CEST 2023
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
function c = geometric_ar1_2d_grid_cell_averaged_cov(lmu,lsd,theta,x,y,dx,dy,varargin)
	% correlation of logarithmic values
	lrfun = @(x,y) exp(-hypot(x,y)/theta);
	% mean of the values
	%mu = exp(lmu+0.5*lsd^2);
	% covariance fucntion of the values
	rfun = @(x,y) logn_corr(lrfun(x,y),lsd,lsd);
	c = cov_cell_averages_2d(rfun,x,y,dx,dy,varargin{:});
end

