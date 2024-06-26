% Wed Jul 11 22:13:14 MSK 2012
% Karl Kästner, Berlin

function [A B p_] = boundary_3d(A, B, Bc, flag)
	% TODO vectorisation and strongly imposing conditions

%	% weakly impose dirichlet boundary conditions
%	for idx=1:size(BC,1)
%	 for jdx=1:size(BC,2)
%		A(BC(idx,jdx),BC(idx,jdx)) = 1e12;
%		B(BC(idx,jdx),BC(idx,jdx)) = 1;
 %        end
%	end

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
end

