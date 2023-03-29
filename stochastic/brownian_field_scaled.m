% Thu 16 Mar 14:21:25 CET 2023
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
%% generate a square (fractal brownian) field where the variance is ellyptic,
%% i.e. increasing at different rates in both axes
%% this is facilitated by cropping and stretching
%% TODO can the kernel directly be adapted?
function e = brownian_field_scaled(H,n,sxy)
		if (length(n)>1)
		if (n(1)~=n(2))
			error('not yet implemented');
		end
		n = n(1);
		end
		
		s2 = sxy.^2;

		scale = max(s2(1)/s2(2),s2(2)/s2(1));
		n_ = ceil(n*scale);
		%n_ = n;

		% on unit square
		e = max(sxy)*brownian_field(0.5,n_);

		% e is a brownian field where the variance increases in all directions at an equal rate
		% interpolate along first direction
		% TODO this should be stochastically interpolated, but the difference is slight, about 1 grid cell
		id = (0:n_-1);
		jd = (0:n-1);
		% if sx(2)/sx(1)=2, and this selects cols 1..n
		e = interp1(id,e,scale*min(1,(s2(1)/s2(2)))*jd,'linear');
		% if sx(2)/sx(1)=2, and this selects rows 1:2:(2*n)
		%jd = (0:(sxy(2)/sxy(1)));
		e = interp1(id,e',scale*min(1,(s2(2)/s2(1)))*jd,'linear')';
		
		% p^2 = sx^2/sy^2 

		% take only the first 1-n
end

