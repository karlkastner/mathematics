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
%
%% function cov_ = geometric_ar1_2d_grid_cell_averaged_cov(lmu,lsd,theta,x,y,dx,dy,varargin)
%%
%% covariance between the grid-cell-averaged values of the continuous ornstein uhlenbeck (ar1) process
function cov_ = ar1_1d_grid_cell_averaged_cov(mu,sd,theta,x,dx,L,varargin)
	% covariance function of the values
	%cfun = @(x) sd*sd*exp(-abs(x)/theta);
	cfun = @(x) sd*sd*(exp(-abs(x)/theta) + exp(-(L-abs(x))/theta));
	% covariance between grid cell averages
	cov_ = cov_cell_averages_1d(cfun,x,dx,varargin{:});
end

