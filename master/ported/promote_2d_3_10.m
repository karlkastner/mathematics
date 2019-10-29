% Tue May 15 17:20:12 MSK 2012
% Karl Kästner
%
% promotes the triangles with tri corner points to triangles with three additional points on the edges
% required to use cubic basis functions
% 
function [P T Bc] = promote_2d_3_10(P, T, Bc)
	% preallocate memory
	lp = size(P,1);
	lt = size(T,1);
	lb = size(Bc,1);
	% each boundary segment is split in three
	Bc = [Bc; zeros(2*lb,3)];
	% there are exactly three new points per triangle and four new points for each boundary
	P = [P; zeros(4*(lt)+lb,2)];

	% there are 7 additional point indices per triangle
	T = [T zeros(size(T,1),7)];

	% triangle egde hashtable
	% there will never be more than lp/2 entries in the hashtable
	import java.util.Hashtable;
	H = java.util.Hashtable(lp);

	% push the boundaries midpoints
	for idx=1:lb
		p1 = min(Bc(idx,1), Bc(idx,2));
		p2 = max(Bc(idx,1), Bc(idx,2));
		% one third points
		P(lp+1,:) = 1/3*(2*P(p1,:) + P(p2,:));
		P(lp+2,:) = 1/3*(P(p1,:) + 2*P(p2,:));
		% split boundary in three
		Bc( idx,:) = [p1   lp+1 Bc(idx,3)];
		Bc(lb+1,:) = [lp+1 lp+2 Bc(idx,3)];
		Bc(lb+2,:) = [lp+2 p2   Bc(idx,3)];
		% memorise the new segment
		H.put(hashkey(p1,p2),lp+1);
		lp = lp+2;
		lb = lb+2;
	end % for idx

	% add edge midpoints to the triangles
	for idx=1:size(T,1)
		p1 = T(idx,1);
		p2 = T(idx,2);
		p3 = T(idx,3);
		% add centre of gravity
		lp = lp+1;
		P(lp,:) = 1/3*(P(p1,:) + P(p2,:) + P(p3,:));
		T(idx,10) = lp;
		% add side points shared with neighbouring triangles or boundaries
		for k1=1:3
			p2_ = min(p2,p3);
			p3_ = max(p2,p3);
			key = hashkey(p2_,p3_);
			val = H.remove(key);
			% this triangle is first and neighbour comes lates
			if (isempty(val))
				% split side into three
				% one third points
				P(lp+1,:) = 1/3*(2*P(p2_,:) + P(p3_,:));
				P(lp+2,:) = 1/3*(P(p2_,:) + 2*P(p3_,:));
				T(idx,3+2*k1-1) = lp+1;
				T(idx,3+2*k1)   = lp+2;
				% memorise new segment
				H.put(key,lp+1);
				lp = lp+2;
			else
				% new side points already created for neighbour triangle
				T(idx,3+2*k1-1) = val(1);
				T(idx,3+2*k1)   = val(1)+1;
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
end % function promote

