% Thu 25 Nov 12:34:20 CET 2021
% Karl Kastner, Berlin
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
%
%% filter input vector with a spatial (two-sided) bandpass in fourier space
%
% function y = bandpass1d_fft(y,fc,p,dx)
function y = bandpass1d_fft(y,fc,p,dx)
	if (isvector(y))
		y = cvec(y);
	end
	n = size(y,1);
	x = dx*(0:n-1)';
	fx = fourier_axis(x);
	normalize = -1;
	S = bandpass1d_discrete_pdf(fx,fc,p,dx,normalize,'f');
	y = real(ifft(((S).^(p/2)).*fft(y)));
end
