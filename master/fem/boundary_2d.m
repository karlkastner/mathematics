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
	p_ = setdiff((1:size(A,1))', p);

	if (nargin() > 3 && 1 == flag)
		% strogly impose bc, remove rows and columns
		A(p,:) = [];
		A(:,p) = [];
		B(p,:) = [];
		B(:,p) = [];
	else	
		% weakly impose the bc
		% todo, without a loop
		for idx=1:length(p)
			% clear boundary rows and columns in the stiffness matrix
			A(p(idx),:) = 0;
			A(:,p(idx)) = 0;
			% set diagonal value close to infinity
			A(p(idx),p(idx)) = 1e12;
			% clear boundary rows and columns in the mass matrix
			B(p(idx),:) = 0;
			B(:,p(idx)) = 0;
			% set diagonal value to 1
			B(p(idx),p(idx)) = 1;
		end
	end
end % boundary_2d

