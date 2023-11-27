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

	if (nargin()<4)
		flag = 0;
	end

	% smoothing
	if (flag)
	S.xy = gaussfilt2(fftshift(Sxy_hat),nf);
	S.xx = gaussfilt2(fftshift(Sxx_hat),nf);
	S.yy = gaussfilt2(fftshift(Syy_hat),nf);
	else
	S.xy = trifilt2(fftshift(Sxy_hat),nf);
	S.xx = trifilt2(fftshift(Sxx_hat),nf);
	S.yy = trifilt2(fftshift(Syy_hat),nf);
	end	

	S.xx = ifftshift(S.xx);
	S.xy = ifftshift(S.xy);
	S.yy = ifftshift(S.yy);

	C = abs(S.xy).^2./(S.xx.*S.yy);
%	C = ifftshift(C);
end

