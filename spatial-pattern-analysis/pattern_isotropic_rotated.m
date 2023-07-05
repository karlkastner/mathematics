% Thu 20 Apr 20:04:29 CEST 2023
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
function [z,x,y] = pattern_isotropic_rotated(n,L,f,af,s,order)

	x = linspace(-L/2,L/2,n)';
	x = repmat(x,1,n);
	y = x';

	% relative coordinates of 4 neighbours
	di = [-0,-1,0,1];
	dj = [-1,0,1, 0];


	% perturbation angle
	ea = 2*pi*rand(n);

	% displacements of grid points
	id = 2:n-1;
	id = 2:n+1;
	dx = zeros(n+2);
	dy = zeros(n+2);

	% neighbour displacements for local rotation
	for idx=1:4
 		dx(id+di(idx),id+dj(idx)) = dx(id+di(idx),id+dj(idx)) + cos(ea)*di(idx)+sin(ea)*dj(idx);
		dy(id+di(idx),id+dj(idx)) = dy(id+di(idx),id+dj(idx)) -sin(ea)*di(idx)+cos(ea)*dj(idx);
	end
	dx_ = dx(2:end-1,2:end-1);
	dx_(:,1) = dx_(:,1) + dx(2:end-1,end);
	dx_(:,end) = dx_(:,end) + dx(2:end-1,1);
	dx_(1,:) = dx_(1,:) + dx(end,2:end-1);
	dx_(end,:) = dx_(end,:) + dx(1,2:end-1);
	dx = dx_;
	dy_ = dy(2:end-1,2:end-1);
	dy_(:,1) = dy_(:,1) + dy(2:end-1,end);
	dy_(:,end) = dy_(:,end) + dy(2:end-1,1);
	dy_(1,:) = dy_(1,:) + dy(end,2:end-1);
	dy_(end,:) = dy_(end,:) + dy(1,2:end-1);
	dy = dy_;

	% smooth the displacements
%	nf = 500;
	for idx=1:4
	dx = lowpass2d_ideal(dx,[1,1],af,order);
	dy = lowpass2d_ideal(dy,[1,1],af,order);
	end
% 	dx = trifilt2(dx,nf);
%	dy = trifilt2(dy,nf);
	
	% rotate coordinates
	x = x+s*dx;
	y = y+s*dy;

	% hexagonal pattern
	z  = 0;
	for idx=0:2
		a = deg2rad(60*idx);
		z = z+cos(2*pi*f*(cos(a)*x + sin(a)*y));
	end
end

