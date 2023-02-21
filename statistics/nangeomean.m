% 2017-06-28 12:57:43.177797184 +0200
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
% geometric mean excluding nan-values
function x = nangeomean(x,varargin)
	x = exp(nanmean(log(x),varargin{:}));
%	if (nargin()<2)
%		dim = 1;
%		if (isvector(x))
%			x = cvec(x);
%		end
%	end
%	n = size(x,dim);
%	s = size(x);
%	s(dim) = 1;
%	y = NaN(s);
%	for idx=1:n
%		if (1==dim)
%			y(1,:) = 
%	end
end

