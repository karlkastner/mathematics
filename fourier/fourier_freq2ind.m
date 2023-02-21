% 2023-02-02 10:06:11.533271532 +0100
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
%
function f = fourier_freq2ind(f,n)
if (0)
	if (mod(n,2) == 0)
	if (max(f(:)) >= n/2)
		max(f(:))
		error('f out of bounds')
	end
	if (min(f(:)) < -n/2)
		min(f(:))
		error('f out of bounds');
	end
	else
	if (max(f(:)) > (n-1)/2)
		max(f(:))
		error('f out of bounds')
	end
	if (min(f(:)) < -(n-1)/2)
		min(f(:))
		error('f out of bounds');
	end
	end
end
	fdx = f < 0;
	f(fdx)  = f(fdx) + n+1;	
	f(~fdx) = f(~fdx) + 1;
%	if (max(f(:))>n || min(f(:)<1))
%		erro('subscript out of range');
%	end
end

