% Mon Apr 30 23:18:31 MSK 2012
% Karl Kästner, Berlin
%
% mark_2d_10
%
% estimates the error of a 2D finite element solution and marks cells for refinement
%
% P : set of points, including points on the boundary, (x, y)
% T : elements (triangles), (p1, p2, p3)
% B : sides coinciding with the domain boundary (p1, p2)
% dV : first order partial derivatives per triangle
% ordering: V must be equally ordered as P, no other assumptions on ordering
% N : iverse mapping from triangle boundaries to triangles
% dV is indexed by triangles

% m : halved order of the pde (1 for second order pdes)
% k : order of basis functions plus 1 (2 for first order basis functions)
% s : order of derivative, which is observed for convergence (0 for solution itself without derivatives)
function [err err_max thresh nH] = mark_2d(mesh, N, dV, m, k, s, mode)
	P   = mesh.P;
	T   = mesh.T;
	Bc  = mesh.Bc;
	area = 0.5*mesh.determinant;
	h_side = mesh.h_side;
	s_angle = mesh.s_angle;
	C = mesh.C;
	err_max = 0;

	% get length
	lt = size(T,1);

	% allocate memory
	err = zeros(lt,1);
	nH = zeros(lt,1);

	% observe convergence of the sth-order derivative
	% see Strang, Fix Theorem 3.7, An Analysis of the Finite Element Method
	rate = min(k-s, 2*(k-m));

	% constants
	C_ = [1/3, 1/13, 1/87, 1/908, 1/12360]; % for Schrödinger mode 1

	switch (mode)
		case {0}
			% C_ = [1/3 1/38 1/870 1/35430 1/882847]; % for poisson mode 0
			C_ = [1/3, 1/13, 1/87, 1/908, 1/12360]; % for Schrödinginger mode 0
			% compute maximum and minimum of the pth-order partial derivatives at the mesh points
			dP_max = zeros(size(P,1),size(dV,2));
			dP_min = zeros(size(P,1),size(dV,2));
			for tdx=1:lt
				dP_max(T(tdx,1),:) = max([dP_max(T(tdx,1),:); dV(tdx,:)]);
				dP_max(T(tdx,2),:) = max([dP_max(T(tdx,2),:); dV(tdx,:)]);
				dP_max(T(tdx,3),:) = max([dP_max(T(tdx,3),:); dV(tdx,:)]);
				dP_min(T(tdx,1),:) = min([dP_min(T(tdx,1),:); dV(tdx,:)]);
				dP_min(T(tdx,2),:) = min([dP_min(T(tdx,2),:); dV(tdx,:)]);
				dP_min(T(tdx,3),:) = min([dP_min(T(tdx,3),:); dV(tdx,:)]);
			end
		case {1}
	end % mode

	% calculate estimated norm of the second derivative per element
	for idx=1:lt
		% seminorm of third derivative
		% for all three triangle neighbours
		switch (mode)
		case {0}
			% take maximum difference over each edge
			nH(idx) = 1/6*sqrt(sum(...
			     ((abs(dP_max(T(idx,2),:) -  dP_min(T(idx,3),:)) + abs(dP_min(T(idx,2),:) -  dP_max(T(idx,3),:)))/h_side(idx,1)).^2 ...
			   + ((abs(dP_max(T(idx,1),:) -  dP_min(T(idx,3),:)) + abs(dP_min(T(idx,1),:) -  dP_max(T(idx,3),:)))/h_side(idx,2)).^2 ...
			   + ((abs(dP_max(T(idx,1),:) -  dP_min(T(idx,2),:)) + abs(dP_min(T(idx,1),:) -  dP_max(T(idx,2),:)))/h_side(idx,3)).^2));
		case {1}
			% Extension to Eriksson, Johnson 1988, Adaptive Finite Element Method for Linear Elliptic Problems
			% fails for s=0 with nonsmooth data
			for jdx=1:3
				b = N(idx,jdx);
				% check that neighbour exists (not a domain boundary)
				if (b > 0)
					dx      = C(idx,1) - C(b,1);
					dy      = C(idx,2) - C(b,2);
					dr_sqr  = dx^2 + dy^2;
					nH_ = 0;
					% take the maximum of the partial derivatives
					for ddx=1:size(dV,2)
						nH_ = max(nH_,abs((dV(idx,ddx) - dV(b,ddx))*dx)/dr_sqr);
						nH_ = max(nH_,abs((dV(idx,ddx) - dV(b,ddx))*dy)/dr_sqr);
					end
					% take maximum maximum over of all three neighbours
					nH(idx) = max(nH(idx), nH_);
				end % if b > 0
			end % jdx
		end % switch mode

		% local error of the element
		err(idx) = C_(k-1) * nH(idx) * (prod(h_side(idx,:))/prod(s_angle(idx,:)))^(1/3*rate);

		% global error over all elements
		err_max = max(err_max, err(idx));
	end % for idx

	% threshold for refinement
	thresh = (0.5^rate)*err_max;
end % err_2d

