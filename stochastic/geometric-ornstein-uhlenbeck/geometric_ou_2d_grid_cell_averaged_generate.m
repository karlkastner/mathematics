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
%% with mean mu_z, standard deviatian sd and correlation length theta_z
%%
%% z     = exp(lz)
%% mean(lz) = mu_z
%% std(lz) = sd_z
%% corr(lz(0,0),lz(x,y)) = exp(-sqrt(x^2 + y^2)/theta_z)
%%
function [z, C, S, lz, lC, lS] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_z,L,n,ni,m,pw)
	% oversampling
	mL = no*L;
	mn = no*n;

	% grid, ordered in fourier style
	[fx,fy] = fourier_axis_2d(mL,mn);
	[x,y] = fourier_axis_2d(mn./mL,mn);
	dx = mL(1)/mn(1);
	dy = mL(2)/mn(2);
	xx = repmat(cvec(x),1,mn(2));
	yy = repmat(rvec(y),mn(1),1);

	% choose sd, theta_z, dx
	% convariance function of the log of the values
	% cov(z), z = log(val)

	% mean of the values, this stays the same, irrespectively of grid cell size
	% mu = exp(mu_z+0.5*sd_z^2)
	mu = lognpdf_mean(mu_z,sd_z);

	% covariance function of the val
	% cfun = @(x,y) logoncorr(lrfun(x,y),mu_z,mu_z,sd_z,sd_z);
	% covariance matrix between grid cell averages
	C = geometric_ou_2d_grid_cell_averaged_cov(mu_z,sd_z,theta_z,xx,yy,dx,dy,ni); %varargin{:});
	S = fft2(C);

	% standard deviation, inverse proportional to grid cell size, due to averaging
	sd = sqrt(C(1,1));

	% correlation matrix
	r = C/C(1,1);

	% the grid cell averages are averages of log-normals and therefore
	% _not_ log-normally distributed, though they are  approximated
	% here by a log normal distribution which matches the first and second
	% moments
	[mu_z_,mu_z__,sd_z_,sd_z__,lC] = lognpdf_moment2par_correlated(mu,mu,sd,sd,r);

	% simulate the log-process
	lz = randn(n);

	% apply window
	if (~isempty(pw))
		w = tukeywin(n(1),pw);
		lCend = lC(1,end/2);
		lC = (w'.*(w.*(lC - lCend))) + lCend;
	end

	% spectral density
	lS = fft2(lC);

	% to oversample the frequency, and still have a periodic pattern within L
	% we have to suppress frequency components with f<1/L

	% TODO set tab n/2 to zero when n is even
	lS = lS(1:o:end,:);
	lS = lS(:,1:o:end);
	S  = S(1:o:end,:);
	S  = S(:,1:o:end);
	lC  = lC(:,[1:L/2,end-L/2+1:end]);
	lC  = lC([1:L/2,end-L/2+1:end],:);
	C  = C(:,[1:L/2,end-L/2+1:end]);
	C  = C([1:L/2,end-L/2+1:end],:);

	% transfer function
	lT = sqrt(lS);

	% correlate (convolve in Fourier domain)
	lz = ifft2(lT.*fft2(lz));

	% remove imaginary part introduced by round off
	lz = real(lz);

	% natural order
	% (note, this is superfluous, as z was iid before convolving, so order does not matter)
	lz = ifftshift(lz);
	
	% make z standard normal, this does not affect the correlation structure
	lz = lz/std(lz(:));

	% scale and translate
	lz = mu_z_ + sd_z_*lz;

	% convert to log-normal
	z = exp(lz);
end % geometric_ou_2d_grid_cell_averaged_generate

