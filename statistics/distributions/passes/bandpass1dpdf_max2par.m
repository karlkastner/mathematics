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
%% function [p] = bandpass1dpdf_mode2par(fc,Sc,p0)
%%
%% transform mode (maxima) of the bandpass spectral density into the paramter
%% of the underlying distribution 
%
function [p] = bandpass1dpdf_max2par(fc,Sc,p0)
	if (nargin()<3)
		p0 = 1;
	end
	% n.b: lsqnonlin works much more reliable than fzero
	p0 = double(p0);
	fc = double(fc);
	Sc = double(Sc);
	p  = lsqnonlin(@(p) bandpass1dpdf_max(fc,abs(p)) - Sc, p0, 0);
end

