% Tue May  8 14:59:52 MSK 2012
% Karl Kästner, Berlin
%
% this applies only homogenous dirichlet boundary conditions at the moment
%
% weakly impose Dirichlet boundary condition
% e.g, u(Gamma) will be 1/1e12 = 10^-12 instead of exactly 0 at the boundary
function [A B p_] = boundary_2d(A,B,Bc,flag)
	% boundary points
	p = unique(Bc(:));
	% non-boundary-points
%	p_ = setdiff((1:size(A,1))', p);
	p_ = (1:size(A,1))';
	p_(p) = [];

	if (nargin() > 3 && 1 == flag)
		% strogly impose bc, remove rows and columns
		A(p,:) = [];
		A(:,p) = [];
		if (nargin() > 1 && ~isempty(B))
			B(p,:) = [];
			B(:,p) = [];
		end
	else	
		% stiffness matrix
		% weakly impose the bc
		% TODO, without a loop
		% clear columns corresponding to boundary points
		for idx=1:length(p)
			A(:,p(idx)) = 0;
		end
		% rows (column operation on transposed matrix is faster)
		A = A';
		for idx=1:length(p)
			A(:,p(idx)) = 0;	
		end
		A = A';
		% set diagonals to a value close to infty
		% TODO scale by 1/h^2
		for idx=1:length(p)
			% set diagonal value close to infinity
			A(p(idx),p(idx)) = 1e12;
		end
	
		% mass matrix
		if (nargin() > 1 && ~isempty(B))
			% clear columns
			for idx=1:length(p)
				B(:,p(idx)) = 0;
			end
			% clear rows
			B = B';
			for idx=1:length(p)
				B(:,p(idx)) = 0;
			end
			B = B';
			% set diagonal value to 1
			for idx=1:length(p)
				B(p(idx),p(idx)) = 1;
			end
		end % if mass matrix supplied
	end
end % boundary_2d

