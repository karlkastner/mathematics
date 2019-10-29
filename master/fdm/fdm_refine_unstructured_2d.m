% Thu May 24 19:13:27 MSK 2012
% Karl Kästner, Berlin

% point coordinates P : [x y boundary_segment], bs = 0, if not on boundary
% neighboors N        : [left right bottom top left2 right2 bottom2 top2]
function [P N] = fdm_refine_unstructured(P, N, Mx, My)
	% preallocate memory
	np = size(P,1);
	P = [P; zeros(np,3)];
	nb = size(B,1);
	B = [B; zeros(nb,1)];

	% refinement in x-direction
	for idx=1:Mx
		a = Mx(idx,1); b = Mx(idx,2);
		% only split, if above and below the old points are no "ghosts"
		if (0 == N(a,7) && 0 == N(a,8) && 0 == N(b, 7) && 0 == N(b, 8)
			% make sure point a is left
			if (P(a,1) > P(b,1))
				c = a;
				a = b;
				b = c;
			end
			% new point
			np = np+1;
			P(np,:) = 0.5*(P(a,:) + P(b,:));
			% link old points to new point
			N( a, 2) = np;
			N( b, 1) = np;
			% link new point to old point
			N(np,1) = a;
			N(np,1) = b;
			% bottom neighbours (regular -> ghost)
			if ()
				% old hanging node exists
				N(np,3) = oh;
				N(oh,4) = np;
				N(oh,8) =  0;
			else
				% no old hanging node, this point becomes hanging node
				N(np,3) = N(a,3);
				N(np,7) = N(b,3);
			end
			% top neighbour
			if ()
				% old hanging node exists
				N(np,4) = oh;
				N(oh,3) = np;
				N(oh,7) =  0;
			else
				% no old hanging node, this point becomes hanging node
				N(np,4) = N(a,4);
				N(np,8) = N(b,4);
			end
	
			% if both points where on the same boundary
			% let this point be on the same boundary, too
			if (P(a,3) > 0 && P(a,3) == P(bm3))
				P(np,3) = P(a,3);
				nb = nb+1;
				B(nb,1) = np;
			end
		else
			% point could not be splitted as at least one of the neighbours needs to be splitted as well
			% save the neighbour(s) and this cell to be splitted (order is important)
			My_(ny_ ) = ...
		end % can be splitted
	end % for idx (Mx)

	% refinement in y-direction
	for idx=1:ny
		a = My(idx,1); b = My(idx,2);
		% only split, if left and right of the old points are no "ghosts"
		if (0 == N(a,5) && 0 == N(a,6) && 0 == N(b, 5) && 0 == N(b, 6)
			% make sure point a is bottom
			if (P(a,2) > P(b,2))
				c = a;
				a = b;
				b = c;
			end
			% new point
			np = np+1;
			P(np,:) = 0.5*(P(a,:) + P(b,:));
			% link old points to new point
			N( a, 4) = np;
			N( b, 3) = np;
			% link new point to old points
			N(np,1) = a;
			N(np,1) = b;

			% left neighbour
			if ()
				% old hanging node exists
				N(np,1) = oh;
				N(oh,2) = np;
				N(oh,6) =  0;
			else
				% no old hanging node, this point becomes hanging node
				N(np,1) = N(a,1);
				N(np,2) = N(b,1);
			end
			% right neighbour
			if ()
				% old hanging node exists
				N(np,2) = oh;
				N(oh,1) = np;
				N(oh,5) =  0;
			else
				% no old hanging node, this point becomes hanging node
				N(np,2) = N(a,2);
				N(np,6) = N(b,2);
			end

			% if both points where on the same boundary
			% let this point be on the same boundary, too
			if (P(a,3) > 0 && P(a,3) == P(bm3))
				P(np,3) = P(a,3);
				nb = nb+1;
				B(nb,1) = np;
			end
		else
			% point could not be splitted as at least one of the neighbours needs to be splitted as well
			% save the neighbour(s) and this cell to be splitted (order is important)
			My_(ny_ ) = ...
			
		end
	end % for idx (My)

	% trim arrays to the right size
	P = P(1:np,:);
	B = B(1:nb,:);
end % fdm_refine_unstructured_2d()

