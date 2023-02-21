% 2021-10-21 09:29:25.938016430 +0200
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
%% integrate the spectral density
%
function I = spectral_density_area(fx,S)
	if (isvector(S))
		S = cvec(S);
	end
	df   = fx(2)-fx(1);
	fdx  = fx>0;
	fdx0 = fx==0;
	if (sum(fdx0)>0)
		I    = (0.5*S(fdx0,:)+sum(S(fdx,:),1))*df;
	else
		I    = sum(S(fdx,:),1)*df;
	end
end

