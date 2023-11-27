% Mon  2 Oct 15:36:18 CEST 2023
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
function [C,S] = coherence_2d(x,y,nf,flag)
	fx  = fft2(x);
	fy  = fft2(y);
	Sxy_hat = fx.*conj(fy);
	Sxx_hat = abs(fx).^2;
	Syy_hat = abs(fy).^2;

	if (nargin()<3)
		nf = 1;
	end

	if (nargin()<4)
		flag = 0;
	end

	% smoothing
%	if (0)
	if (nf > 0)
		S.xy = ifftshift(gaussfilt2(fftshift(Sxy_hat),nf));
		S.xx = ifftshift(gaussfilt2(fftshift(Sxx_hat),nf));
		S.yy = ifftshift(gaussfilt2(fftshift(Syy_hat),nf));
	end
%	else
%	S.xy = trifilt2(fftshift(Sxy_hat),nf);
%	S.xx = trifilt2(fftshift(Sxx_hat),nf);
%	S.yy = trifilt2(fftshift(Syy_hat),nf);
%	end
%	end	

	% TODO rename into radial average
	Sr.xx = periodogram_radial(S.xx);
	Sr.xy = periodogram_radial(S.xy);
	Sr.yy = periodogram_radial(S.yy);

	C = abs(Sr.xy.mu).^2./(Sr.xx.mu.*Sr.yy.mu);
%	C = ifftshift(C);
end

