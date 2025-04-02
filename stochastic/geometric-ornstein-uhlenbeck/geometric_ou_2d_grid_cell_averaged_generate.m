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
%% simulate a geometric two dimensional spatial Ornsten-Uhlenbeck process,
%% i.e. the process with circular boundary conditions:
%%
%% z     = exp(lz)
%% mean(lz) = mu_z
%% std(lz) = sd_z
%% corr(lz(0,0),lz(x,y)) = exp(-sqrt(x^2 + y^2)/theta_z)
%%
%% where z of the continuos process is log normally distributed
%% with mean mu_z, standard deviatian sd_z and correlation length theta_z
%%
%% the discrete process is only close to the continuous process for intermediate
%% correlation length where dx << theta << L
%%
%% artefacts for not small dx/theta are reduced by grid-cell averaging (oversampling in the spatial domain)
%% artefacts for not large L/theta are reduced by oversampling in the spectral domain and windowing
%%
%% when approximated in the midpoint scheme (m_spatial = 1), i.e. without oversampling
%% in space, the OU process is approximated by an two-dimensional spatial AR1 process
%%
%% input:
%%	mu_z    : mean of the continuous and discrete process
%%	sd_z    : standard deviation of the continuous process,
%%		  approached by the discrete process in the limit dx -> 0
%%	theta_z : correlation length of the continuous process,
%%		  approached by the discrete process in the limit df = 1/L -> 0
%%
%% L            : spatial extent
%% n = L/dx     : number of grid points
%% 
%% m_spatial    : number of oversampling points along one dimension when oversampling the spectral density
%% m_spectral   : number of oversampling points along one dimension when oversampling the correlation
%%
%% pwC          : tukey-window parameter in the spatial domain
%%		  spectral window not necessary, as C is sampled without osciallation
function [z, C, S, out, lz, lC, lS, w] = geometric_ou_2d_grid_cell_averaged_generate(mu_z,sd_z,theta_z,L,n,m_spatial,m_spectral,pwC)
	% oversampling in spectral domain
	mL = m_spectral*L;
	mn = m_spectral*n;

	% grid, ordered in fourier style
	[fx,fy,frr] = fourier_axis_2d(mL,mn);
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
	[mu_lz, sd_lz] = lognpdf_moment2par(mu_z,sd_z);
	% for dx/theta->0
	theta_lz = geometric_ou_correlation_length_lin2log(sd_lz,theta_z);

	% covariance function of the val
	% cfun = @(x,y) logoncorr(lrfun(x,y),mu_z,mu_z,sd_z,sd_z);
	% covariance matrix between grid cell averages
	C = geometric_ou_2d_grid_cell_averaged_cov(mu_lz,sd_lz,theta_lz,xx,yy,dx,dy,m_spatial);

	% apply window
	if (~isempty(pwC) && (pwC>0))
		w = tukeywin_circular(mn,pwC);
		Cend = C(1,end/2);
		C = (w.*(C - Cend)) + Cend;
	else
		w = 0;
	end
	S = fft2(C);

	% standard deviation of the discrete process
	% inverse proportional to grid cell size, due to averaging
	sd_zd = sqrt(C(1,1));

	% correlation matrix of the discrete process
	r_zd = C/C(1,1);

	% the grid cell averages are averages of log-normals and therefore
	% _not_ log-normally distributed, though they are  approximated
	% here by a log normal distribution which matches the first and second
	% moments
	[mu_lzd,mu_lzd,sd_lzd,sd_lzd,lC] = lognpdf_moment2par_correlated(mu_z,mu_z,sd_zd,sd_zd,r_zd);

	% simulate the log-process
	lz = randn(n);
	% either we have to subtract the mean here or st lS(0) to zero,
	% as the mean can be scaled considerably when theta is large
	lz = lz - mean(lz,'all');

	% spectral density
	lS = fft2(lC);

	% to oversample the frequency, and still have a periodic pattern within L
	% we have to suppress frequency components with f<1/L
	lS = resample_density_by_averaging_2d(lS,m_spectral);
	S  = resample_density_by_averaging_2d(S,m_spectral);

	if (m_spectral>1)
		% the covariance structure is recomputed here, as it is modified
		% by oversampling in the spectral domain
		lC = ifft2(lS);
		C  = ifft2(S);
	end	

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
	lz = mu_lzd + sd_lzd*lz;

	% convert to log-normal
	z = exp(lz);

	out.mu_lz = mu_lz;
	out.sd_lz = sd_lz;
	out.theta_lz = theta_lz;
end % geometric_ou_2d_grid_cell_averaged_generate

