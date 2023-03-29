% Sat 11 Jun 17:31:52 CEST 2022
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
function [b,xy,S,f,R] = anisotropic_pattern(L,n,f0,sxy,mode)
	if (nargin()<5)
		mode = 's1d';
	end

	S  = struct();

	% coordinate axes in real space
	xy.x = linspace(0,L(1),n(1));
	xy.y = linspace(0,L(2),n(2));
	% coordinate axes in frequency_space
	f.x  = fourier_axis(xy.x);
	f.y  = fourier_axis(xy.x);

	switch (lower(mode))
	case {'exact'}
		%if (L(1) ~= L(2) || n(1) ~= n(2))
		%	error('not yet implemented');
		%end

		e = brownian_field_scaled(0.5,n,sxy)';

		if (nargout()>2)
			[S.S2d,xy.x,xy.y,R] = anisotropic_pattern_pdf(L,n,f0,sxy);
			[f.x,f.y] = fourier_axis_2d(L,n);
		end

		
		% scale from 1x1 to L*L square
		e = sqrt(L(1))*e;

		% scale it by stdev
		%e = sxy(1)*e;

		phase = e + rvec(xy.x)*f0;

		b = cos(2*pi*phase);
	case {'s1d'}
		% approximate, via transfer function
		S.x = phase_drift_pdf(f.x,f0,sxy(1),true);
		S.y = phase_drift_parallel_pdf(f.y,sxy(2));
		S.S2d = cvec(S.y)*rvec(S.x);
		if (nargout()>4)
			R = ifft(S.S2d);
			R = R/R(1);	
		end

		% white noise to be filtered into pattern
		e  = randn(n);
		fe = fft2(e);

		% transfer function
		T = sqrt(S.S2d);

		% fourier transform of pattern
		fb = T.*fe;
	
		% periodogram of pattern
		S.Shat = abs(fb).^2;

		% pattern	
		b = ifft2(fb);

		% this should be real, but due to roundoff error, have a small complex part
		b = real(b);
	case {'s2d'}
		[S.S2d,xy.x,xy.y] = anisotropic_pattern_pdf(L,n,f0,sxy);
		if (nargout()>4)
			[R] = anisotropic_pattern_acf(L,n,f0,sxy);
		end
	
		e  = randn(n);
		fe = fft2(e);
		T  = sqrt(S.S2d);
		fb = T.*fe;
		b  = ifft2(fb);
		b  = real(b);
	otherwise 
		error('here')
	end
end

