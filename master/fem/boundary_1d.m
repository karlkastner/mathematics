% Tue Jun 19 14:33:10 MSK 2012
% Karl Kästner, Berlin

function [A B p_] = boundary_1d(A,B,Bc,bcflag)
	% boundary points
	p = Bc(:,1);
	% non-boundary-points
	p_ = setdiff((1:size(A,1))', p);

	if (nargin() > 3 && 1 == bcflag)
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
			% set diagonla value to 1
			B(p(idx),p(idx)) = 1;
		end
	end
end % end

