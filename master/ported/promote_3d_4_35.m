% Thu Jul 12 15:47:38 MSK 2012
% Karl Kästner, Berlin

%  1  2  3  4  5  6  7
%  1  3  6 10 15 21 28
%  1  4 10 20 35
%  0  0  0
%  0  0  0 

% the element has 20 points, 1 point more than necessary for the sake of symmetry
function [P T Bc] = promote_3d_4_35(P, T, Bc)
	lt1 = size(T,1);
	lb1 = size(Bc,1);
	np1 = size(P,1);

	% allocate memory
	Bc = [Bc zeros(lb1,12)]; 
	T  = [T  zeros(lt1,31)];
	P  = [P; zeros(31*lt1,3)];
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
			if (isempty(val))
				val = np1;
				H.put(key,val);
				P(val+1,:) = 1/4*(3*P(p1,:) +   P(p2,:));
				P(val+2,:) = 1/4*(2*P(p1,:) + 2*P(p2,:));
				P(val+3,:) = 1/4*(  P(p1,:) + 3*P(p2,:));
				np1 = np1+3;
			end
			Bc(idx,3+3*jdx-2) = val+1;
			Bc(idx,3+3*jdx-1) = val+2;
			Bc(idx,3+3*jdx  ) = val+3;
		end
		% boundary face inner points
		p1 = Bc(idx,1);
		p2 = Bc(idx,2);
		p3 = Bc(idx,3);
		% this key cannot be in the hash yet
		key = hashkey(p1,p2,p3);
		H.put(key,np1);
		P(np1+1,:) = 0.125*(6*P(p1,:) +   P(p2,:) +   P(p3,:));
		P(np1+2,:) = 0.125*(  P(p1,:) + 6*P(p2,:) +   P(p3,:));
		P(np1+3,:) = 0.125*(  P(p1,:) +   P(p2,:) + 6*P(p3,:));
		Bc(idx,13) = np1+1;
		Bc(idx,14) = np1+2;
		Bc(idx,15) = np1+3;
		np1 = np1+3;
	end
	
	% split the inner edges into two
	% local side indices
	q = nchoosek(1:4,2);
	% local face indices
	f = nchoosek(1:4,3);
	for idx=1:lt1
		% edge points
		for jdx=1:6
			p1 = T(idx,q(jdx,1));
			p2 = T(idx,q(jdx,2));
			key = hashkey(p1,p2);
			val = H.get(key);
			if (isempty(val))
				val = np1;
				H.put(key, val);
				P(val+1,:) = 1/4*(3*P(p1,:) +   P(p2,:));
				P(val+2,:) = 1/2*(  P(p1,:) +   P(p2,:));
				P(val+3,:) = 1/4*(  P(p1,:) + 3*P(p2,:));
				np1 = np1+3;
			end
			T(idx,4+3*jdx-2) = val+1;
			T(idx,4+3*jdx-1) = val+2;
			T(idx,4+3*jdx  ) = val+3;
		end
		% face inner points
		for kdx=1:4
			p1 = T(idx,f(kdx,1));
			p2 = T(idx,f(kdx,2));
			p3 = T(idx,f(kdx,3));
			key = hashkey(p1,p2,p3);
			val = H.get(key);
			if (isempty(val))
				val = np1;
				H.put(key,val);
				P(val+1,:) = 0.125*(6*P(p1,:) +   P(p2,:) +   P(p3,:));
				P(val+2,:) = 0.125*(  P(p1,:) + 6*P(p2,:) +   P(p3,:));
				P(val+3,:) = 0.125*(  P(p1,:) +   P(p2,:) + 6*P(p3,:));
				np1 = np1+3;
			end
			T(idx,22+3*kdx-2) = val+1;
			T(idx,22+3*kdx-1) = val+2;
			T(idx,22+3*kdx  ) = val+3;
		end
		% element mid point
		p1 = T(idx,1);
		p2 = T(idx,2);
		p3 = T(idx,3);
		p4 = T(idx,4);
		np1 = np1+1;
		P(np1,:) = 0.25*( P(p1,:) + P(p2,:) + P(p3,:) + P(p4,:) );
		T(idx,35) = np1;
	end % for idx
	% other than in 2D edges are shared by more than 2 tetras and a
	% simple consistency check is not possible here
	P = P(1:np1,:);
end % function promote_3d_4_20

