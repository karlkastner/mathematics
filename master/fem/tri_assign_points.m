function pdx = tri_assign_points(mesh,X,Y)
	% points inside the triangle form a convex combination
	T = mesh.T;
	P = mesh.P;
	for idx=1:size(T,1)
		% set up vandermonde matrix
		A = [             1,             1,             1;
                      P(T(idx,1),1), P(T(idx,2),1), P(T(idx,3),1);
                      P(T(idx,1),2), P(T(idx,2),2), P(T(idx,3),2)];
		% interpolation constant for each point
		c = A \ [ones(1,length(X)); X(:)'; Y(:)' ];
		% points with 0 <= c(:) <= 1 are convex combinations
		tmp = [0 <= c(1,:) & 0 <= c(2,:) & 0 <= c(3,:) ...
                                c(1,:) <= 1 & c(2,:) <= 1 & c(3,:) <= 1];	
		pdx{idx} = find(tmp);
	end % for idx
end % function tri_assign_points

