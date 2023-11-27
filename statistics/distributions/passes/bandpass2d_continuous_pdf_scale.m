% Fri 22 Apr 16:00:51 CEST 2022
% Karl Kästner, Berlin
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
% scale (normalization) factor of the spectral density of the 2D bandpass filter
% the scale factor is equal to the maximum of the density
% mode (maximum) of the sd of the 2d bp
% inverse of the normalization constant
function [Sc] = bandpass2d_continuous_pdf_mode(L,n,a,order)
	fr = fourier_axis(L,n);
	fr = fr(fr>=0);
	S  = bandpass2d_continuous_pdf(fr,a,order);
	df  = 1/L;
	ciS = sum(S)*df;
	Sc  = 1./ciS;
end

