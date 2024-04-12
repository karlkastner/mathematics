% Wed 21 Dec 11:27:03 CET 2022
% Karl Kästner, Berlin
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
%% smooth (filter) the  2D image z with a gaussian window
%% apply periodic boundary conditions
% function [zbar, dof, w] = gaussfilt2(z, nf)
function [zbar, dof, w] = gaussfilt2(z, nf)
	
	[w, dof] = gausswin2(size(z), nf);

	% Fourier transform of window
	fw = fft2(w);

	% by definition, the ft is real
	fw = real(fw);

	% apply mean filter by convolution theorem
	zbar = ifft2(fw.*fft2(z));

	if (isreal(z))
		zbar = real(zbar);
	end
end

