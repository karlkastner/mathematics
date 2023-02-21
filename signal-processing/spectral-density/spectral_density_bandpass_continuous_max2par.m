% Fri  7 Jan 16:11:57 CET 2022
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
%% transform mode (maxima) of the bandpass spectral density into the paramter
%% of the underlying distribution 
%
% function [p] = spectral_density_bandpass_max2par(fc,Sc,p0)
function [p] = spectral_density_bandpass_continuous_max2par(fc,Sc,p0,pp)
	if (nargin()<3)
		p0 = 1;
	end
	if (nargin()<4)
		pp = [];
	end
%	p  = fzero(@(p) spectral_density_bandpass_continuous_scale(fc,abs(p),pp) - Sc, p0);
	% n.b: lsqnonlin works much more reliable than fzero
	p  = lsqnonlin(@(p) spectral_density_bandpass_continuous_scale(fc,abs(p),pp) - Sc, p0);
%	p  = fzero(@(p) log(spectral_density_bandpass_continuous_scale(fc,abs(p),pp)) - log(Sc), p0);
	p  = abs(p);
end

