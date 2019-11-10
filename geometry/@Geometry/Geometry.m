% Mi 11. Nov 10:25:11 CET 2015
% Karl Kastner, Berlin

% TODO split this class -> own class for triangles
% TODO own class for polygons
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

classdef Geometry
	properties
	end % properties
	methods (Static)
		S                = arclength(x,y,z);
		[xb, s, t]       = base_point(x0,x1,x2);
		[xb, s, t]       = base_point_limited(x0,x1,x2);
		[cx, cy]         = centroid(x,y);
		[c,R]            = curvature(varargin);
		d2               = distance2(xy1,xy2);
		d                = distance(xy1,xy2);
		cosa             = dot(dx1,dy1,dx2,dy2);
		[cosa, g, H]     = ddot(a,b,c);
		H                = edge_length(E,P);
		[cosa, alpha]    = enclosed_angle(a,b);
		[X, Y, elem, s]  = enclosing_triangle(X,Y);
		xy	         = hexagon(xy0,scale,alpha);
		[in]		 = inTetra(X,Y,Z,X0,Y0,Z0);
		[in, flag]       = inTetra2(P1,P2,P3,P4,P0); 
		in               = inPolygon(xp,yp,x,y);
		[px, py]         = intersect(x1,y1,x2,y2);
		[flag, c]        = inTriangle(X,Y,X0,Y0);
		cflag 		 = quad_isconvex(X,Y);
		[flag, s, t, p, q, den] = lineintersect(p1,p2,q1,q2);
		[x0, y0]         = mittenpunkt(x,y);
		[x0, y0]         = nagelpoint(x,y);
		[flag, c]        = onLine(X,x0);
		[x0, y0]         = orthocentre(x,y);
		area             = poly_area(x,y);
                [val, flag]      = poly_set(poly,pval,x0,y0,val);
		E                = poly_edges(xp,yp);
		W 	         = poly_width(X,Y);
		[X0, Y0, FDX, DEN]  = polyxpoly(X,Y);
		[xy0p,id,p,dis]        = project_to_curve(xyc,xy0);
		[xp, yp, p, d]   = plumb_line(x1,y1,x2,y2,x0,y0,cflag);
		[P0, Rc, R]      = random_simplex(varargin);
		v                = sphere_volume(r);
		V                = tetra_volume(X,Y,Z);
		c                = tobarycentric(X,Y,X0,Y0)
		c                = tobarycentric1(P1,P2,P0);
		c                = tobarycentric2(P1,P2,P3,P0);
		c                = tobarycentric3(P1,P2,P3,P4,P0);
		[xc, yc]         = tri_centroid(x,y);
		cosa 		 = tri_angle(x,y);
		[A]              = tri_area(X,Y);
		s                = tri_semiperimeter(X,Y);
		[h2, h]		 = tri_height(X,Y);
		[h2, h]		 = tri_edge_length(X,Y);
		[x, y]           = tri_edge_midpoint(X,Y);
		[d2, d]          = tri_distance_opposit_midpoint(X,Y);
		[Xc, Yc, R]      = tri_excircle(X,Y);
%		h	         = tri_height(X,Y);
		[x0, y0, R]      = tri_incircle(X,Y);
		[isacute, Xc, Yc] = tri_isacute(X,Y);
		[isobtuse, Xc, Yc] = tri_isobtuse(X,Y);
		[d2, d] 	 = tri_side_length(X,Y)
	end % methods (Static)
end % Geometry

