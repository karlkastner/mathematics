% Fr 13. Nov 18:02:32 CET 2015
% Karl Kastner, Berlin
%% centroidal voronoi tesselation
% TODO improve triangulation by inserting points at the centre of short edges (element lmin/lmax)
function [X Y mesh] = cvt(xb,yb,n,abstol)
	xmax = max(xb);
	xmin = min(xb);
	ymax = max(yb);
	ymin = min(yb);
	% 1) sample randomly initial points in desired area
	% 1a) sample random points
	X = xmin+(xmax-xmin)*rand(n,1);
	Y = xmin+(ymax-ymin)*rand(n,1);
	% 1b) remove points outside of the domain
	[in] = inpolygon(X,Y,xb,yb);
	X = X(in);
	Y = Y(in);


	n = length(X);
	jdx=1;
	% while convergence criterium not met
	while (1)
		% convex hull
		X(end+1:end+4) = [-1e3 1e3 1e3 -1e3];
		Y(end+1:end+4) = [-1e3 -1e3 1e3 1e3];

		% 2) create voronoi polygon
		[V,C] = voronoin([X Y]);
%		fdx = isinf(V) & V > 0;
%		V(fdx) = 1e3;
%		fdx = isinf(V) & V < 0;
%		V(fdx) = -1e3;

%		figure(jdx);
%		clf
%		for idx=1:n
%			c = C{idx};
%			c(end+1) = c(1);
			%patch(V(c,1),V(c,2),idx);
%			hold on
%			plot(V(c,1),V(c,2))
%		end
%		axis([xmin,xmax,ymin,ymax]);
%		plot(X,Y,'.');
		
		% (- cut polygons to region of interest)
		% 3) compute centre of mass for each polygon
		% 4) move points to centre of mass
		% discard last four points
		X(end-3:end) = [];
		Y(end-3:end) = [];

		Xold = X;
		Yold = Y;
		% except the convex hull points
		for idx=1:n
			c = C{idx};
			c(end+1) = c(1);
			[X(idx), Y(idx)] = poly_centroid(V(c,1),V(c,2));
		end
		% overshoot
	%	X = Xold + 2*(X - Xold);
	%	Y = Yold + 2*(Y - Yold);

		% 5) project points to domain boundary or constraint position on the domain
		X = max(xmin,min(xmax,X));
		Y = max(ymin,min(ymax,Y));
%		plot(X,Y,'.');
		if ( max(hypot(X-Xold,Y-Yold)) <= abstol)
			break;
		end
		jdx=jdx+1;
	end
	printf('Converged in %d iterations\n',jdx-1);
	% generate mesh
	l = max(cellfun(@length,C));
	elem = NaN(length(C),l);
	for idx=1:length(C)
		c = C{idx};
		elem(idx,1:length(c)) = c;
	end
	mesh = Mesh();
	% quick fix
	V = min(V,1e3);
	V = max(V,-1e3);
	mesh.point = V;
	mesh.elem = elem;
end


