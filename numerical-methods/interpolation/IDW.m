% Sun Jul  6 11:02:33 WIB 2014
% Karl Kastner, Berlin

%% spatial averaging by inverse distance weighting
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
classdef IDW < Interpolator
	properties
		% no additional properties
	end % properties

	methods
	function obj = IDW(Rmax,Rmin,order,nverify,aspect_ratio)
		% call superclass constructor
		obj = obj@Interpolator(Rmax,Rmin,order,nverify,aspect_ratio);
	end

	function [vt et ot fdx res dx dy obj] = interpolate_(obj,Xs,Vs,x0,Rmin,cfunc)
		% find all points within radius Ri
		% TODO use a quadtree here
		% TODO is bsxfunc faster?
		D = obj.dist(Xs,repmat(x0,size(Xs,1),1));
		D2 = (D(:,1).^2 + D(:,2).^2);
		R2max = obj.Rmax*obj.Rmax;
		fdx = find(D2 < R2max & D2 > Rmin*Rmin);
		%dx = D(fdx,1);
		%dy = D(fdx,2);
		if (~isempty(fdx))
			% unnormalised weights
			w = 1./sqrt(D2(fdx));
			% normalise weights
			w = w/sum(w);
			vt  = w'*Vs(fdx,:);
			% The error estimate is not found in the literature
			% and somewhat experimental
			et = sqrt(w'*((Vs(fdx,:) - vt).^2));
			% actually the order is not well defined
			ot = 0;
		else
			vt = NaN;
			et = NaN;
			ot = NaN;
		end
		res = [];
		dx = [];
		dy = [];	
		% constant interpolaton
		%if (d2min <= R2max)
		%else
		%	
		%end
	end % interpolate_
	
	end % methods
end % classdef
