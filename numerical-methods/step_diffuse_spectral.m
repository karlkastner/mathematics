% Mon  2 May 14:18:38 CEST 2022
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% analytic solution to the heat equation
%
% dy/dt = e*D^2*y
%
% fft    : (analytic) solution via fourier transform
% no fft : analytic solution by convolution
% note   : for non-smooth solutions, variant the no-fft assures positivity
%          the "no-fft" variant is also implimented via fourier transform
% function zt = diffuse_spectral(t,z0,dx,e)
function zt = step_diffuse_spectral(t,z0,n,dx,e)
	% fundamental solution:
	% 1/(4*pi*t)^(n/2)*exp(-r.^2/(4*t))

%	if (nargin()>3&&~isempty(e))
%		t = e*t;
%	end
	e = e(1);

	if (nargin()<5)
		usefft = true;
	end

	ndim = length(n);
%	if (isvector(z0))
%		ndim = 1;
%	else
%		ndim = 2;
%	end

	% convolve with Gaussian in fourier space
	switch (ndim)
	case {1}
		%n = length(z0);
		L_ = n.*dx;
		or = 2*pi*fourier_axis(L_,n);
		or = cvec(or);
		fz0 = fft(z0);
		fg = exp(-(or.*or)*e*t);
%clf
%h%old on
		%plot(ifft(fg.^100))
%		plot((fg.^100),'--')
%pause
		fz = fg.*fz0;
		zt = ifft(fz);
	case {2}
		%n = size(z0);
		L_ = n./dx;
		[fx, fy, fr, ft] = fourier_axis_2d(L_,n);
		or = 2*pi*fr;
		%dx = rvec(L)./rvec(n);
		fz = fft2(z0);
		fg = exp(-(or.*or).*e*t);
		fz = fg.*fz;
		zt = ifft2(fz);
		%s = sum(zt(:));
		%zt = zt./s;
	otherwise
		error('not yet implemented');
	end
end % diffuse_spectral

