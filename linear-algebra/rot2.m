% Fr 15. Jan 09:55:03 CET 2016
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
%% rotation matrix from angle
% function R = rot2(alpha,idx,jdx,n)
function R = rot2(alpha_rad,idx,jdx,n)
	if (nargin()<2)
		idx=1;
	end
	if (nargin()<3)
		jdx=2;
	end
	if (nargin() < 3)
		n=max(idx,jdx);
	end
%	idx = min(idx_,jdx_);
%	jdx = max(idx_,jdx_);
if (~issym(alpha_rad))
	R = speye(n);
else
	R = sym(eye(n));
end
	c = cos(alpha_rad);
	s = sin(alpha_rad);
	R(idx,idx) =  c;
	R(jdx,jdx) =  c;
	R(idx,jdx) = -s; % signs were swapped
	R(jdx,idx) = +s;
end

