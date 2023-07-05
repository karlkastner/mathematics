% Sun 11 Jul 21:39:38 CEST 2021
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
% function y = lowpass1d_fft(x,rho,order)
function y = lowpass1d_fft(x,rho,order)
	if (nargin()<3)
		order = 1;
	end
	% order_fft = sqrt(order_implicit)
	order = sqrt(order);
	if (isvector(x))
		x = cvec(x);
	end
	n  = size(x,1);
	fx = fourier_axis(1,n);
	dx = 1/n;
	%S  = spectral_density_lowpass_discrete(fx,rho,order,dx);
	S  = lowpass1d_discrete_pdf(fx,rho,order,dx);
	y  = ifft(S.^(order/2).*fft(x));

end

