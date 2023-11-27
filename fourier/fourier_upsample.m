% Wed  8 Nov 09:26:56 CET 2023
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
function x = fourier_upsample(x,nu)
	if (isvector(x))
		x = cvec(x);
	end
	n = size(x);
	f = fft(x);
	% alternatively, instead of removing the imaginary n/2+1 component,
	% we can insert its complex conjugate
	if (mod(n(1),2) == 0)
	f = [f(1:n(1)/2,:);
             zeros(nu-n(1)+1,n(2));
	     f(n(1)/2+2:end,:)];
	else
	n(1) = n(1)+1;
	f = [f(1:n(1)/2,:);
             zeros(nu-n(1)+1,n(2));
	     f(n(1)/2+1:end,:)];
	end
	f = nu/n(1)*f;
	x = ifft(f);	
end

