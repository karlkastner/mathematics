% Mon May  7 22:04:00 MSK 2012
% Karl Kästner
%
% promotes the triangles with tri corner points to triangles with three additional points on the edges
% required to use quadratic basis functions
% 
function [P T Bc] = promote_2d_3_6(P, T, Bc)
	% preallocate memory
	lp = size(P,1);
	lt = size(T,1);
	lb = size(Bc,1);
	% each boundary segment is split in two
	Bc = [Bc; zeros(lb,3)];
	% there are exactly one new point for each triangle and boundary segment less by one
	%P = [P; zeros(lt+lb,2)];
	P = [P; zeros(2*lt+lb+1-lp,2)];

	% there are at most twice as many points plus one
	T = [T zeros(size(T))];
	
	% triangle egde hashtable
	% there will never be more than lp/2 entries in the hashtable
	import java.util.Hashtable;
	H = java.util.Hashtable(lp);

	% push the boundaries midpoints
	for idx=1:lb
		p1 = Bc(idx,1);
		p2 = Bc(idx,2);
		% midpoint
		lp = lp+1;
		P(lp,:) = 0.5*(P(p1,:) + P(p2,:));
		H.put(hashkey(p1,p2),lp);
		% split boundary in two
		Bc(idx,2) = lp;
		lb = lb+1;
		Bc(lb,:) = [lp p2 Bc(idx,3)];
	end % for idx

	% add edge midpoints to the triangles
	for idx=1:size(T,1)
		p1 = T(idx,1);
		p2 = T(idx,2);
		p3 = T(idx,3);
		for k1=1:3
			key = hashkey(p2,p3);
			val = H.remove(key);
			% this triangle is first and neighbour comes lates
			if (isempty(val))
				lp = lp+1;
				P(lp,:) = 0.5*(P(p2,:)+P(p3,:));
				H.put(key,lp);
				T(idx,3+k1) = lp;
			else
				% point already created for neighbour triangle
				T(idx,3+k1) = val(1);
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

