% Thu Jul 12 15:47:38 MSK 2012
% Karl Kästner, Berlin


% the element has 20 points, 1 point more than necessary for the sake of symmetry
function [P T Bc] = promote_3d_4_20(P, T, Bc)
	lt1 = size(T,1);
	lb1 = size(Bc,1);
	np1 = size(P,1);

	% allocate memory
	Bc = [Bc zeros(lb1, 6)]; 
	T  = [T  zeros(lt1,16)];
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
%				['honk' num2str(val)]
				Bc(idx,3+2*jdx-1) = val+1;
				Bc(idx,3+2*jdx  ) = val+2;
			else
				P(np1+1,:) = 1/3*(2*P(p1,:) +   P(p2,:));
				P(np1+2,:) = 1/3*(  P(p1,:) + 2*P(p2,:));
				H.put(key,np1);
				Bc(idx,3+2*jdx-1) = np1+1;
				Bc(idx,3+2*jdx  ) = np1+2;
%				['m' num2str(val)]
				np1 = np1+2;
			end
		end
		% boundary face midpoint
		p1 = Bc(idx,1);
		p2 = Bc(idx,2);
		p3 = Bc(idx,3);
		np1 = np1+1;
		P(np1,:) = 1/3*(P(p1,:) + P(p2,:) + P(p3,:));
		key = hashkey(p1,p2,p3);
		% key is cannot be in hash yet
		H.put(key,np1);
		Bc(idx,10) = np1;
	end
	
	% split the inner edges into two
	% local side indices
	q = nchoosek(1:4,2);
	% local face indices
	f = nchoosek(1:4,3);
	for idx=1:lt1
		% process all sides
		for jdx=1:6
			p1 = T(idx,q(jdx,1));
			p2 = T(idx,q(jdx,2));
			key = hashkey(p1,p2);
			val = H.get(key);
			if (~isempty(val))
				T(idx,4+2*jdx-1) = val+1;
				T(idx,4+2*jdx  ) = val+2;
			else
				P(np1+1,:) = 1/3*(2*P(p1,:) +   P(p2,:));
				P(np1+2,:) = 1/3*(  P(p1,:) + 2*P(p2,:));
				H.put(key, np1);
				T(idx,4+2*jdx-1) = np1+1;
				T(idx,4+2*jdx  ) = np1+2;
				np1 = np1+2;
			end
		end
		% midpoints of the faces
		for jdx=1:4
			p1 = T(idx,f(jdx,1));
			p2 = T(idx,f(jdx,2));
			p3 = T(idx,f(jdx,3));
			key = hashkey(p1,p2,p3);
			val = H.get(key);
			if (isempty(val))
				np1 = np1+1;
				P(np1,:) = 1/3*(P(p1,:) + P(p2,:) + P(p3,:));
				H.put(key,np1);
				T(idx,16+jdx) = np1;
			else
				T(idx,16+jdx) = val;
			end
		end
	end
	% other than in 2D edges are shared by more than 2 tetras and a
	% simple consistency check is not possible here
end % function promote_3d_4_20

