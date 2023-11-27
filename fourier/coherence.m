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
function C = mycoherence(x,y,nf)
	fx  = fft(x);
	fy  = fft(y);
	Sxy_hat = fx.*conj(fy);
	Sxx_hat = abs(fx).^2;
	Syy_hat = abs(fy).^2;

	% smotthing
	Sxy = gaussfilt1(Sxy_hat,nf);
	Sxx = gaussfilt1(Sxx_hat,nf);
	Syy = gaussfilt1(Syy_hat,nf);
	
	C = abs(Sxy).^2./(Sxx.*Syy);
end

