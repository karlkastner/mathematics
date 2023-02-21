% Wed  9 Aug 11:56:50 CEST 2017
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
%% weighted root mean square
% varargin can be dimension
%
%     mu = 1/(sum w) sum w x
%      e = x - mu
%
% function wrms = wrms(w,x,varargin)
function wrms = wrms(w,x,varargin)
	if (isvector(x))
		x = cvec(x);
	end
	if (isvector(w))
		w = cvec(w);
	end
	wms  = sum(bsxfun(@times,w,x.^2),varargin{:})./sum(w,varargin{:});
	wrms = sqrt(wms);
end
