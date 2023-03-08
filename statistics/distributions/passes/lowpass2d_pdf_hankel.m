% Sat  4 Mar 13:44:15 CET 2023
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
%
%% spectral density of the two-dimensional lowpass filter with autocorrelation
%%
%% r = exp(-a*sqrt(x^2 + y^2))
%%
%% efficiently estimated with gauss-laguerre integration and 1D-FFT:
%%
%% for a radially symmetric function, the radial density is
%% Sr(r) = S2d(r,0) = S2d(0,r)
%% with density S2d and autocorrelation R2d
%% S2d = F_2d^-1 (R2d)
%% by the slicing theorem:
%% S2d(r,0) = F_1d^-1 (int R2d(x,y) dy)
function [Sr,Rbar] = lowpass2d_pdf_hankel(fr,a,order,m);
	if (nargin()<4)
		% this is about 4-digits arrurate
		m = 20;
	end
 	[w,fy] = int_1d_gauss_laguerre(m);
	fy    = fy/a;
	% average aucorrelation
	% the scaling with 2 pi a of the integral can be omitted,
	% as the density is scaled/normalized later
	Rbar = exp(-a*hypot(cvec(fr),rvec(fy)))*(exp(a*cvec(fy)).*w);
	% radial density
	Sr = abs(real(fft(Rbar)));
	% scale to 1
	Sr = Sr/Sr(1);
	if (nargin()>2 && ~isempty(order))
		Sr = Sr.^order;
	end
	% TODO normalize
end

