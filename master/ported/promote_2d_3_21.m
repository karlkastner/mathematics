% Sun May 20 22:21:10 MSK 2012
% Karl Kästner
%
% promotes the triangles with tri corner points to triangles with three additional points on the edges
% required to use quadratic basis functions
% 
%		TODO interior point location not optimal - see java port where it is already corrected
%
function [P T Bc] = promote_2d_3_21(P, T, Bc)
	% preallocate memory
	lp = size(P,1);
	lt = size(T,1);
	lb = size(Bc,1);
	% each boundary segment is split into five
	Bc = [Bc; zeros(4*lb,3)];
	% there are exactly 10 new points per triangle and four new points for each boundary
	P = [P; zeros(12*lt+2*lb,2)];
	% there are 18 additional point indices per triangle
	T = [T zeros(size(T,1),18)];

	% triangle egde hashtable
	% there will never be more than lp/2 entries in the hashtable
	import java.util.Hashtable;
	H = java.util.Hashtable(lp);

	% split boundaries into five parts with four new points each
	for idx=1:lb
		p1 = min(Bc(idx,1), Bc(idx,2));
		p2 = max(Bc(idx,1), Bc(idx,2));
		% 1/5th points
		P(lp+1,:) = 0.2*(4*P(p1,:) +   P(p2,:));
		P(lp+2,:) = 0.2*(3*P(p1,:) + 2*P(p2,:));
		P(lp+3,:) = 0.2*(2*P(p1,:) + 3*P(p2,:));
		P(lp+4,:) = 0.2*(  P(p1,:) + 4*P(p2,:));
		% split boundary in four
		Bc( idx,:) = [p1   lp+1 Bc(idx,3)];
		Bc(lb+1,:) = [lp+1 lp+2 Bc(idx,3)];
		Bc(lb+2,:) = [lp+2 lp+3 Bc(idx,3)];
		Bc(lb+3,:) = [lp+3 lp+4 Bc(idx,3)];
		Bc(lb+4,:) = [lp+4 p2   Bc(idx,3)];
		% memorise the new segment
		H.put(hashkey(p1,p2),lp);
		lp = lp+4;
		lb = lb+4;
	end % for idx

	% promote each triangle
	for idx=1:size(T,1)
		p1 = T(idx,1);
		p2 = T(idx,2);
		p3 = T(idx,3);
		% add 6 interior points
		P(lp+1,:) = 0.10*(8*P(p1,:) +   P(p2,:) +   P(p3,:));
		P(lp+2,:) = 0.10*(  P(p1,:) + 8*P(p2,:) +   P(p3,:));
		P(lp+3,:) = 0.10*(  P(p1,:) +   P(p2,:) + 8*P(p3,:));
		P(lp+4,:) = 0.05*(2*P(p1,:) + 9*P(p2,:) + 9*P(p3,:));
		P(lp+5,:) = 0.05*(9*P(p1,:) + 2*P(p2,:) + 9*P(p3,:));
		P(lp+6,:) = 0.05*(9*P(p1,:) + 9*P(p2,:) + 2*P(p3,:));
		T(idx,16:21) = [lp+(1:6)];
		lp = lp+6;
		% add side points shared with neighbouring triangles or boundaries
		for k1=1:3
			p2_ = min(p2,p3);
			p3_ = max(p2,p3);
			key = hashkey(p2_,p3_);
			val = H.remove(key);
			% this triangle is first and neighbour comes later
			if (isempty(val))
				% split side into five
				% 1/5th points
				P(lp+1,:) = 0.2*(4*P(p2_,:) +   P(p3_,:));
				P(lp+2,:) = 0.2*(3*P(p2_,:) + 2*P(p3_,:));
				P(lp+3,:) = 0.2*(2*P(p2_,:) + 3*P(p3_,:));
				P(lp+4,:) = 0.2*(  P(p2_,:) + 4*P(p3_,:));
				T(idx,3+4*k1-3) = lp+1;
				T(idx,3+4*k1-2) = lp+2;
				T(idx,3+4*k1-1) = lp+3;
				T(idx,3+4*k1  ) = lp+4;
				% memorise new segment
				H.put(key,lp);
				lp = lp+4;
			else
				% new side points already created for neighbour triangle
				T(idx,3+4*k1-3) = val(1)+1;
				T(idx,3+4*k1-2) = val(1)+2;
				T(idx,3+4*k1-1) = val(1)+3;
				T(idx,3+4*k1  ) = val(1)+4;
			end
			% circular shift
			tmp = p1;
			p1 = p2;
			p2 = p3;
			p3 = tmp;
		end
	end % for idx

	% verify edge consistency
	if (~H.isEmpty())
		K = H.keySet().toArray();
		for idx=1:length(K)
			[(K(idx) / 2^23) mod(K(idx), (2^23-1)) H.get(K(idx))];
		end
		error('promote_2d','Inconsistent Mesh');
	end
end % function promote_2d_3_21

