% Thu May 24 19:28:38 MSK 2012
% Karl Kästner, Berlin

% fist derivative
% set up the 5-point convection operator for an unstructured FDM 2D grid
% central differences

% point coordinates P : [x y boundary_segment], bs = 0, if not on boundary
% neighboors N        : [left right bottom top left2 right2 bottom2 top2]
% c : constant coefficient [cx cy]
function A = fdm_assemble_d1_2d(P, N, c)
	np = size(P,1);
	buf = zeros(9*np,3);
	nb = 0;
	% for each point
	for idx=1:np 
		% top, bottom and left,right can be swapped without harm
		% x steps
		hx = P(N(idx,1),1) - P(N(idx,2),1);
		% y steps
		hy = P(N(idx,3),2) - P(N(idx,4),2);
	
		% main diagonal is zero
	
		% x-entries
	
		% left
		if (0 == N(idx,5))
			% left neighbour exists
			nb = nb+1;
			buf(nb,:)   = [idx, N(idx,1), c(1)/hx];
		else
			% left neighbour is interpolated
			buf(nb+1,:) = [idx, N(idx,1), 0.5*c(1)/hx];
			buf(nb+2,:) = [idx, N(idx,5), 0.5*c(1)/hx];
			nb = nb+2;
		end
		% right
		if (0 == N(idx,6))
			nb = nb+1;
			buf(nb,:)   = [idx, N(idx,2), c(1)/hx];
		else
			buf(nb+1,:) = [idx, N(idx,1), 0.5*c(1)/hx];
			buf(nb+2,:) = [idx, N(idx,6), 0.5*c(1)/hx];
			nb = nb+2;
		end
	
		% y-entries
	
		% bottom
		if (0 == N(idx,7))
			nb = nb+1;
			buf(nb,:) = [idx, N(idx,3), c(2)/hy];
		else
			buf(nb+1,:) = [idx, N(idx,3), c(2)/hy];
			buf(nb+2,:) = [idx, N(idx,7), c(2)/hy];
			nb = nb+2;
		end
	
		% top
		if (0 == N(idx,8))
			nb = nb+1;
			buf(nb,:)   = [idx, N(idx,4), c(2)/hy];
		else
			buf(nb+1,:) = [idx, N(idx,4), 0.5*c(2)/hy];
			buf(nb+2,:) = [idx, N(idx,8), 0.5*c(2)/hy];
			nb = nb+2;
		end
	end % for idx
	
	A = sparse(buf(1:nb,1), buf(1:nb,2), buf(1:nb,2));
end % fdm_assemble_d1_2d()

% a minimal initial grid for the unstructured FDM

	P = [ kron(linspace(0, L0(1), nx)', ones(ny,1)), ...
              kron(ones(nx,1), linspace(0, L0(2), ny)') ];
	N = zeros(nx*ny,9);
	% set left neighbours
	% set right neighbours
	% set 


