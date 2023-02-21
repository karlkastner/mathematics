% 2015-04-07 18:35:59.479232653 +0200
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
%
%% mid point between neighbouring vector elements
%
function x = mid(x,dim)
	if (isvector(x))
		x = 0.5*(x(1:end-1)+x(2:end));
	else
		if (nargin()<2)
			dim = 1;
		end
		switch (ndims(x))
		case {2}
			if (1 == dim)
				x = 0.5*(x(1:end-1,:)+x(2:end,:));
			else
				x = 0.5*(x(:,1:end-1)+x(:,2:end));
			end
		case {3}
			switch (dim)
			case {1}
				x = 0.5*(x(1:end-1,:,:)+x(2:end,:,:));
			case {2}
				x = 0.5*(x(:,1:end-1,:)+x(:,2:end,:));
			case {3}
				x = 0.5*(x(:,:,1:end-1)+x(:,:,2:end));
			end
		case {4}
			switch (dim)
			case {1}
				x = 0.5*(x(1:end-1,:,:,:)+x(2:end,:,:,:));
			case {2}
				x = 0.5*(x(:,1:end-1,:,:)+x(:,2:end,:,:));
			case {3}
				x = 0.5*(x(:,:,1:end-1,:)+x(:,:,2:end,:));
			case {4}
				x = 0.5*(x(:,:,:,1:end-1)+x(:,:,:,2:end));
			end
			
		end
	end
end
