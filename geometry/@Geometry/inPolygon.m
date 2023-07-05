% Fr 27. Nov 17:42:09 CET 2015
% Karl Kastner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% flag points contained in a polygon
%% much faster than matlab internal function
% function in = inpolygon(xp,yp,x,y)
function in = inpolygon(xp,yp,x,y)
	x  = cvec(x);
	y  = cvec(y);
	xp = cvec(xp);
	yp = cvec(yp);

	% slice domain vertically
	xmax = max(xp);
	xmin = min(xp);
	ymin = min(yp);
	Lx   = xmax-xmin;

	% size of vertical slice
	np = length(xp);
	dx = Lx/np;

	E = Geometry.poly_edges(xp,yp);

	% index of each polygon point to its respective vertical
	sp = floor((xp-xmin)/dx)+1;
	% indices for each edge
	mE = [min([sp(E(:,1)) sp(E(:,2))],[],2), max([sp(E(:,1)) sp(E(:,2))],[],2)];
	sE = arrayfun(@(a,b) colon(a,b),mE(:,1),mE(:,2),'uniformoutput',false);

	spi = cell(np+1,1);
	% inverse index from colum to points, this is not unique
	% TODO this is kind of expensive
	for idx=1:length(sE)
		for jdx=1:length(sE{idx})
			spi{sE{idx}(jdx)}(end+1) = idx;
		end
	end

%	for idx=1:length(spi)
%		figure(3)
%		clf
%		q1 = [xp(E(spi{idx},1)) yp(E(spi{idx},1))]';
%		q2 = [xp(E(spi{idx},2)) yp(E(spi{idx},2))]';
%		plot([q1(1,:); q2(1,:)], [q1(2,:); q2(2,:)],'.-');
%		pause
%	end
	
	% determine slices of input points
	sdx = floor((x-xmin)/dx)+1;

	% for each point in test set
	c = zeros(length(x),1);
	for idx = 1:length(x)
		% get vertical up
		vdx = sdx(idx);
		if (vdx > 0 && vdx < length(spi))
		
		% get polygon points in the vertical
		% count how often the vertical north of this point cuts the vertical
		p1 = [x(idx);y(idx)];
		p2 = [x(idx);ymin];
		q1 = [xp(E(spi{vdx},1)) yp(E(spi{vdx},1))]';
		q2 = [xp(E(spi{vdx},2)) yp(E(spi{vdx},2))]';
		[flag s t] = Geometry.lineintersect(p1,p2,q1,q2);
%		figure(3)
%		clf
%		plot(xp,yp,'.k')
%		hold on
%		plot([q1(1,:); q2(1,:)], [q1(2,:); q2(2,:)],'.-');
%		plot([p1(1,1),p2(1,1)],[p1(2,1),p2(2,1)],'.-');
%		pause
	
		% to avoid problems to count exact hits on points twice or not at all,
		% this has to be gt and lte
		% this is still a problem with zero length edges,
		% the assumption is that all edges have a non-zero length
		c(idx) = sum(s > 0 & t > 0 & t <= 1);
		end
	end
	% points which cut the boundary an even number of times are inside the polygon
	in = logical(mod(c,2));
end

