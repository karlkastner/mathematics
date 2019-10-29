% Sun Aug 12 13:58:32 MSK 2012
% Karl Kästner, Berlin

%Multigrid matrices
%- in promote move BC after T, save for each point the tetrahedron it originates from
%- P -> T_mesh -> T_tree -> parent in T_tree-1 -> T_mesh-1 (problem, not unique) (index)
%
%promotion routine:
%	C_j : test function matrix of T index in other mesh for point P_i from original mesh
%	P_i' * C U = phi_interp U
%	=> phi_interp is a row vector, put into row i, columns T(adx)

function A = mg_mat(P,T,Q,poly)
	tol = 1e-7;
	lt2 = size(T,2);
	np = size(P,1);
	nq = size(Q,1);

	for idx=1:nq;
		% load coordinates of point to be interpolated
		q = [1 Q(idx,:)]';
		% fetch index of the triangle containing the point in the current mesh
		tdx = M(idx);
		% fetch points of this triangle
		A  = P(T(tdx,:),:);
		Va = vander_3d(A,poly);
		% check that q is contained in the tetra (convex combination, all coefficients positive)
		c = inv(A(1:3,1:3))*q;
		if (abs(sum(c-abs(c))) > tol)
			error('mg_mat','wrong index');
		end
		% get test functions
		C = inv(A);
		% get interpolation coeffcients
		% q' C U = u(p)
		phi = q'*C;
		for adx=1:lt2
			% avoid to put zero values into the sparse matrix
			if (abs(phi(adx) > tol))
				buf(nb,1) = idx;
				buf(nb,2) = T(tdx,adx);
				buf(nb,3) = phi(adx);
				nb=nb+1;
			end
		end
	end

	A = spdiags(buf(:,1), buf(:,2), buf(:,3), nq, np);
end

