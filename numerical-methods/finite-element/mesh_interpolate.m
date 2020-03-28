% Fri Mar 14 13:49:10 WIB 2014
% Karl Kastner, Berlin

% uses the finite elment method to interpolate the values V, X, Y to the corner
% points mesh.P

function mesh = mesh_interpolate(mesh, X, Y, V, pfile)
	javaaddpath('/home/pia/Documents/master/thesis/src/fem/class');
	javaaddpath('/usr/share/java/jama.jar');

	% translate to avoid round off errors
	xm = mean(mesh.P(:,1));
	ym = mean(mesh.P(:,2));
	mesh.P(:,1) = mesh.P(:,1) - xm;
	mesh.P(:,2) = mesh.P(:,2) - ym;
	X = X - xm;
	Y = Y - ym;

	switch(mesh.type)
		case {'triangle'}
			if (nargin() < 5)
				pdx_  = tri_assign_points(mesh,X,Y);
				% convert to matlab array
				for idx=1:length(pdx_)
					pdx{idx} = arrayfun(@(x) x, pdx_(idx).toArray);
				end
				save('pdx.mat','pdx');
			else
				load('pdx.mat','pdx');
			end
			mesh = tri_interpolate(mesh,pdx,X,Y,V);
		case {'quadrilateral'}
			pdx  = quad_assign_points(mesh,X,Y);
			mesh = quad_interpolate(mesh,pdx,X,Y,V);
		otherwise
			error('unsupported mesh type');
	end
	% undo translation (todo store the values and do not invert +-)
	mesh.P(:,1) = mesh.P(:,1) + xm;
	mesh.P(:,2) = mesh.P(:,2) + ym;
end % interpolate

function pdx = tri_assign_points(mesh,X,Y)
	% if object creation fails and the class path exists,
	% than verify that the object was compiled with the same java version
	% and recompile with "javac -source 1.6 -target 1.6 Mesh_2d.java"
	Bc = [0 0];
	javamesh = javaObject('Mesh_2d', mesh.P, mesh.T, zeros(0,2) )
	pdx = javamesh.assign_points(X,Y);
end

% constant interpolation (per element)
function mesh = tri_interpolate(mesh,pdx,X,Y,V)
	% TODO, linear interpolation
%	mesh = quad_interpolate(mesh,pdx,X,Y,V);
	% set up the discretisation matrix
%	for idx=1:size(mesh.T,1)
%		% construct the vandermonde matrix
%		A = [1 1 1;
%			P(T(idx,1),1) P(T(idx,2),1) P(T(idx,),1)
%			P(T(idx,1),2) P(T(idx,2),2) P(T(idx,),2) ];
%		% regression matrix
%		A = [ ones(length(pdx{idx}),1) X(pdx{idx}) Y(pdx({idx}) ];
%		b = V(pdx{idx});
%		% interpolation coefficients
%		c = A \ b;
%		c = A \ b;
%	end

	% set up interpolation matrix
	A = interpolation_matrix(mesh,pdx,X,Y);
	% remove points, which are not defined
	s = sum(A>0,2);
	fdx = find(s > 2);
	% strip rows of ill defined points
	A = A(fdx,:);
	Vi = NaN(size(mesh.P,1),1);
	Vi(fdx) = (A') \ V;
	E_ = (A')*Vi(fdx) - V;
	for idx=1:length(pdx)
		E(idx) = norm(E_(pdx{idx}));
	end
	mesh.V = Vi;
	mesh.E = E;
end % tri_interpolate_points

function pdx = quad_assign_points(mesh,X,Y)
	P = mesh.P;
	T = mesh.T;
	% assign points to the quadrilaterals
	% for each quadrilateral, find points in the bounding box
	for idx=1:size(T,1)
		xmax = max(P(T(idx,:),1));
		xmin = min(P(T(idx,:),1));
		ymax = max(P(T(idx,:),2));
		ymin = min(P(T(idx,:),2));
		pdx{idx} = find( X >= xmin & X <= xmax  & Y >= ymin & Y <= ymax);
		% TODO, check that points are really within the quadrilateral
	end % for idx
%		% TODO : this is a stupip bounding box approach, better use interpolation coefficients
%		l = min(x_p);
%		r = max(x_p);
%		b = min(y_p);
%		t = max(y_p);
%		fdx = find(x >= l && x <= r && y >= b && y <= t);

end % quad_assign_points

function mesh = quad_interpolate(mesh,pdx,X,Y,V)

	% set up interpolation matrix
	A = interpolation_matrix(mesh,pdx,X,Y);
	% remove points, which are not defined
	s = sum(A>0,2);
	fdx = find(s > 2);
	% strip rows of ill defined points
	A = A(fdx,:);
	Vi = NaN(size(mesh.P,1),size(V,2));
	Vi(fdx,:) = (A') \ double(V);
	E_ = (A')*Vi(fdx,:) - V;
	for idx=1:length(pdx)
		E(idx) = norm(E_(pdx{idx}))/sqrt(length(pdx{idx}));
	end
	mesh.V = Vi;
	mesh.E = E;

	% constant interpolation
	% TODO linear interpolation
%	V = [];
%	E = [];
%figure(3);
%clf
%	m = size(V_,2);
%	for idx=1:size(mesh.T,1)
%		vi = V_(pdx{idx},:);
%		if (0 == length(vi))
%			V(idx,1:m) = NaN;
%			%E(idx,1:m) = NaN;
%		else
%			% interpolated depth / altitude
%			V(idx,1:m) = mean(vi,1);
%			% interpolation error
%			E(idx,1:m) = 0; %norm(vi - V(idx,1))/sqrt(length(vi));
%		end % if
%		px = mesh.P(:,1);
%		py = mesh.P(:,2);
%		plot(X(A),Y(A),'k.')
%		hold on
%	end % for idx
%	mesh.V = V;
%	mesh.E = E;
end % function quad_interpolate

