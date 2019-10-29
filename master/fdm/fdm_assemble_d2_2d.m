% Thu May 24 19:15:52 MSK 2012
% Karl Kästner, Berlin

% second derivative - diffusion
% set up the 5-point Laplacian operator for an unstructured FDM 2D grid
% point coordinates P : [x y boundary_segment], bs = 0, if not on boundary
% neighboors N        : [left right bottom top left2 right2 bottom2 top2]
% c : constant coefficient [cx cy]
function A = fdm_assemble_d2_2d(P, N, c)
	np = size(P,1);
	buf = zeros(9*np,3);
	nb = 0;
	% for each point
	for idx=1:np 
		% top, bottom and left,right can be swapped without harm
		% x steps
		hxl = P(N(idx,1),1) - P(idx,1);
		hxr = P(N(idx,2),1) - P(idx,1);
		hxc = hxl+hxr;
		% y steps
		hyl = P(N(idx,3),2) - P(idx,2);
		hyr = P(N(idx,4),2) - P(idx,2);
		hyc = hxl+hxr;
	
		% main diagonal
		nb = nb+1;
		buf = -2*c(1)/(hxl*hxr) - 2*c(2)/(hyl*hyr)
	
		% x-entries
	
		% left
		if (0 == N(idx,5))
			% left neighbour exists
			nb = nb+1;
			buf(nb,:)   = [idx, N(idx,1), 2*c(1)/(hxl*hxc)];
		else
			% left neighbour is interpolated
			buf(nb+1,:) = [idx, N(idx,1), c(1)/(hxl*hxc)];
			buf(nb+2,:) = [idx, N(idx,5), c(1)/(hxl*hxc)];
			nb = nb+2;
		end
		% right
		if (0 == N(idx,6))
			nb = nb+1;
			buf(nb,:)   = [idx, N(idx,2), 2*c(1)/(hxr*hxc)];
		else
			buf(nb+1,:) = [idx, N(idx,1), c(1)/(hxr*hxc)];
			buf(nb+2,:) = [idx, N(idx,6), c(1)/(hxr*hxc)];
			nb = nb+2;
		end
	
		% y-entries
	
		% bottom
		if (0 == N(idx,7))
			nb = nb+1;
			buf(nb,:)   = [idx, N(idx,3), 2*c(2)/(hyl*hyc)];
		else
			buf(nb+1,:) = [idx, N(idx,3), c(2)/(hyl*hyc)];
			buf(nb+2,:) = [idx, N(idx,7), c(2)/(hyl*hyc)];
			nb = nb+2;
		end
	
		% top
		if (0 == N(idx,8))
			nb = nb+1;
			buf(nb,:)   = [idx, N(idx,4), 2*c(2)/(hyr*hyc)];
		else
			buf(nb+1,:) = [idx, N(idx,4), c(2)/(hyr*hyc)];
			buf(nb+2,:) = [idx, N(idx,8), c(2)/(hyr*hyc)];
			nb = nb+2;
		end
	end % for idx
	
	A = sparse(buf(1:nb,1), buf(1:nb,2), buf(1:nb,2));
end % fdm_assemble_d2_2d

