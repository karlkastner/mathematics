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
%% simulate a grid cell averaged stochastic process exp(z)
%% where z follows a geometric Ornstein-Uhlenbeck (AR1) process
%% with mean lmu, standard deviatian sd and stationary autocorrelation
%%
%% val = exp(z)
%% mean(z) = lmu
%% std(z)  ) ls
%% corr(z(0,0),z(x,y)) = exp(-sqrt(x^2 + y^2)/theta)
%%
%%
function val = geometric_ar1_2d_grid_cell_averaged_generate(lmu,lsd,theta,L,n,varargin)
	% grid, ordered in fourier style
	[x,y] = fourier_axis_2d([1,1],n);
	dx = L(1)/n(1);
	dy = L(2)/n(2);
	x = x*dx;
	y = y*dy;
	xx = repmat(cvec(x),1,length(y));
	yy = repmat(rvec(y),length(x),1);

	% choose sd, theta, dx
	% convariance function of the log of the values
	% cov(z), z = log(val)

	% mean of the values, this stays the same, irrespectively of grid size
	% mu = exp(lmu+0.5*lsd^2)
	mu = logn_mean(lmu,lsd);

	% covariance function of the values
	%cfun = @(x,y) logn_corr(lrfun(x,y),lmu,lmu,lsd,lsd);
	% covariance matrix between grid cell averages
	cov_ = geometric_ar1_2d_grid_cell_averaged_cov(lmu,lsd,theta,xx,yy,dx,dy,varargin{:});

	% standard deviation, decreasing with grid cell size, due to averaging
	sd = sqrt(cov_(1,1));

	% correlation matrix
	r = cov_/cov_(1,1);

	% the grid cell averages are averages of log-normals and therefore
	% _not_ log-normally distributed, though they are  approximated
	% here by a log normal distribution which matches the first and second
	% moments
	[lmu_,lsd_,lr_] = logn_moment2par_correlated(mu,sd,r);
	% simulate the process
	z = randn(n);

	% correlate (convolve in fourier domain)
	z = ifft2(fft2(lr_).*fft2(z));
	% scale and translate
	val = exp(lmu_ + lsd_*z);
	% natural order
	val = ifftshift(val);
end

