% Thu Jul 12 15:47:38 MSK 2012
% Karl Kästner, Berlin

function [P T Bc] = promote_3d_4_10(P, T, Bc)
	lt1 = size(T,1);
	lb1 = size(Bc,1);
	np1 = size(P,1);

	% allocate memory
	Bc = [Bc zeros(lb1,3)]; 
	T  = [T zeros(lt1,6)];
	% TODO, preallocate points

	import java.util.Hashtable;
	H = java.util.Hashtable(6*lt1);

	% split edges of the boundary segments 
	% TODO, save identification number of the boundary
	% local side indices
	q = nchoosek(1:3,2);
	for idx=1:lb1
		% process local sides
		for jdx=1:3
			p1 = Bc(idx,q(jdx,1));
			p2 = Bc(idx,q(jdx,2));
			key = hashkey(p1,p2);
			val = H.get(key);
			if (~isempty(val))
				Bc(idx,3+jdx) = val;
			else
				np1 = np1+1;
				P(np1,:) = 0.5*(P(p1,:) + P(p2,:));			
				H.put(key,np1);
				Bc(idx,3+jdx) = np1;
			end
		end
	end
	
	% split the inner edges into two
	% local side indices
	q = nchoosek(1:4,2);
	for idx=1:lt1
		% process all sides
		for jdx=1:6
			p1 = T(idx,q(jdx,1));
			p2 = T(idx,q(jdx,2));
			key = hashkey(p1,p2);
			val = H.get(key);
			if (~isempty(val))
				T(idx,4+jdx) = val;
			else
				np1 = np1+1;
				P(np1,:) = 0.5*(P(p1,:) + P(p2,:));			
				H.put(key,np1);
				T(idx,4+jdx) = np1;
			end
		end
	end
	% other than in 2D edges are shared by more than 2 tetras and a
	% simple consistency check is not possible here
end % function promote_3d_4_10

