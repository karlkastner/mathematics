% Wed 29 Mar 16:30:11 CEST 2023
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
%% simulate a grid cell averaged stationary stochastic process exp(z)
%% where z follows a geometric Ornstein-Uhlenbeck (AR1) process
%% with mean lmu, standard deviatian sd and correlation length theta
%%
%% val     = exp(z)
%% mean(z) = lmu
%% std(z)) = lsd
%% corr(z(0,0),z(x,y)) = exp(-sqrt(x^2 + y^2)/theta)
%%
function [val,cov_] = ar1_1d_grid_cell_averaged_generate(mu,sd,theta,L,n,varargin)
	% grid, ordered in fft style
	x = fourier_axis(n./L,n);
	dx = L(1)/n(1);

	% choose sd, theta, dx
	% convariance function of the log of the values
	% cov(z), z = log(val)

	% mean of the values, this stays the same, irrespectively of grid size
	% mu = exp(lmu+0.5*lsd^2)

	% covariance matrix between grid cell averages
	% this is effectively just the reshaped first row of the full covariance
	% matrix, which has identical rows (up to shift)
	% note, the second argument must be indeed inf otherwise theta is doubled
	cov_ = ar1_1d_grid_cell_averaged_cov(mu,sd,theta,x,dx,inf,varargin{:});

	% standard deviation, decreasing with grid cell size, due to averaging
	sd = sqrt(cov_(1));

	% correlation matrix
	r = cov_/cov_(1);

	% simulate the process
	z = randn(n,1);

	% correlate (convolve in Fourier domain)
	z = ifft(fft(r).*fft(z));

	% natural order 
	% (note, this is superfluous, as z was iid before convolving, so order does not matter)
	% compared to conv, this is even incorrect
	% as z is periodic,cutting it in half and flipping the halves does not change anything
	z = fftshift(z);
	
	% make z standard normal, this does not affect the correlation structure
	z = z/std(z(:));

	% scale and translate
	val = mu + sd*z;

end % ar1_1d_grid_cell_averaged_generate


