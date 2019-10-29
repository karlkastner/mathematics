% Wed Jun 20 23:01:14 MSK 2012
% Karl Kästner, Berlin

%classdef Refine_2d
%   properties
%  		P
%		lp
%		T
%		N
%		lt
%		Bc
%		lb
%		L
%		Lx
%   end % properties
%   methods
%      function obj = Refine_2d(obj,P,T,N,Bc)
%		obj.P = P;
%		obj.T = T;
%		obj.N = N;
%		obj.Bc = Bc;
%		obj.lp = size(P,1);
%		obj.lt = size(T,1);
%		obj.lb = size(Bc,1);
%      end

% we skip the obj-stuff in matlab here to avoid the unpleasant "obj.member"

 function [P T Bc N void] = refine_2d_21(P, T, Bc, N, M, varargin)
	void = [];
	lp = size(P,1);
	lt = size(T,1);
	lb = size(Bc,1);

	% compute side length
	L = sqrt( [(P(T(:,2),1) - P(T(:,3),1)).^2 + (P(T(:,2),2) - P(T(:,3),2)).^2, ...
	           (P(T(:,1),1) - P(T(:,3),1)).^2 + (P(T(:,1),2) - P(T(:,3),2)).^2, ...
	           (P(T(:,1),1) - P(T(:,2),1)).^2 + (P(T(:,1),2) - P(T(:,2),2)).^2 ]);

	% find longest sides
	% Attention: this leads to non optimal results if the
	% longest side is not uinque
	Lx = 1 +  (L(:,2) > L(:,1)).*(L(:,2) > L(:,3)) ...
		+ 2*(L(:,3) > L(:,1)).*(L(:,3) > L(:,2));

	% push boundaries
	Sh = java.util.Hashtable(lt);
	for bdx=1:length(Bc)
		Sh.put(hashkey(Bc(bdx,1), Bc(bdx,2)), [-bdx 0]);
	end
	
	% push sides of old triangles here
	for tdx=1:lt
		key = hashkey(T(tdx,1), T(tdx,2));
		val = Sh.put(key,[tdx 3])';
		if (~isempty(val))
			Sh.put(key,[val tdx 3]);
		end
		key = hashkey(T(tdx,1), T(tdx,3));
		val = Sh.put(key,[tdx 2])';
		if (~isempty(val))
			Sh.put(key,[val tdx 2]);
		end
		key = hashkey(T(tdx,2), T(tdx,3));
		val = Sh.put(key,[tdx 1])';
		if (~isempty(val))
			Sh.put(key,[val tdx 1]);
		end
	end

	Ph = java.util.Hashtable(lp);



	af = zeros(lt,1);
	% for each triangle to be refined
	for mdx=1:length(M)
		[P lp Bc lb Ph Sh af] = split_edges(M(mdx),P,lp,T,N,Bc,lb,Ph,Sh,L,Lx,af);
	end

%		keySet = Sh.keySet();
%		K = keySet.toArray();
%		for kdx = 1:length(K)
%			key = K(kdx);
%			% remove hanging node from hash
%			val = Sh.get(key)'
%		end

	% for each triangle
	% can be simplified by noting the affected elements
	F = find(af);
%	for tdx=1:lt
	for fdx=1:length(F)
		tdx=F(fdx);
		%[T lt] = partition(tdx,T,lt,N,Ph);
		[T lt N Sh] = partition(tdx,T,lt,N,Ph,L,Lx,Sh);
	end

end % refine_2d

function [P lp Bc lb Ph Sh af] = split_edges(tdx,P,lp,T,N,Bc,lb,Ph,Sh,L,Lx,af)
	af(tdx) = 1;
	p1 = T(tdx,1);
	p2 = T(tdx,2);
	p3 = T(tdx,3);
	% remove edges
	[T Sh] = remove_edges(tdx,T,Sh);
	% split its three sides
	key = hashkey(p2,p3);
	val = Ph.get(key);
	if (isempty(val))
		lp = lp+1;
		P(lp,:) = 0.5*(P(p2,:) + P(p3,:));
		Ph.put(key,lp);
		% recursively split the longest side of the neighbouring triangle
		if (N(tdx,1) > 0)
			[P lp Bc lb Ph Sh af] = split_longest(N(tdx,1),P,lp,T,N,Bc,lb,Ph,Sh,L,Lx,af);
		else
			[Bc lb Sh] = split_boundary(-N(tdx,1),Bc,lb,lp,Sh);
		end
	end
	key = hashkey(p1,p3);
	val = Ph.get(key);
	if (isempty(val))
		lp = lp+1;
		P(lp,:) = 0.5*(P(p1,:) + P(p3,:));
		Ph.put(key,lp);
		if (N(tdx,2) > 0)
			[P lp Bc lb Ph Sh af] = split_longest(N(tdx,2),P,lp,T,N,Bc,lb,Ph,Sh,L,Lx,af);
		else
			[Bc lb Sh] = split_boundary(-N(tdx,2),Bc,lb,lp,Sh);
		end
	end
	key = hashkey(p1,p2);
	val = Ph.get(key);
	if (isempty(val))
		lp = lp+1;
		P(lp,:) = 0.5*(P(p1,:) + P(p2,:));
		Ph.put(key,lp);
		if (N(tdx,3) > 0)
			[P lp Bc lb Ph Sh af] = split_longest(N(tdx,3),P,lp,T,N,Bc,lb,Ph,Sh,L,Lx,af);
		else
			[Bc lb Sh] = split_boundary(-N(tdx,3),Bc,lb,lp,Sh);
		end
	end
end % split_edges()

function [P lp Bc lb Ph Sh af] = split_longest(tdx,P,lp,T,N,Bc,lb,Ph,Sh,L,Lx,af)
	af(tdx) = 1;
	[T Sh] = remove_edges(tdx,T,Sh);
	% get longest side
	switch(Lx(tdx))
		case {1}
			p1 = T(tdx,2);
			p2 = T(tdx,3);
			n3 = N(tdx,1);
		case {2}
			p1 = T(tdx,3);
			p2 = T(tdx,1);
			n3 = N(tdx,2);
		case {3}
			p1 = T(tdx,1);
			p2 = T(tdx,2);
			n3 = N(tdx,3);
	end
	key = hashkey(p1,p2);
	val = Ph.get(key);
	% longest side was not yet split
	if (isempty(val))
		lp = lp+1;
		P(lp,:) = 0.5*(P(p1,:) + P(p2,:));
		Ph.put(key,lp);
		if (n3 > 0)
			% recursive splitting of longest edge
			[P lp Bc lb Ph Sh af] = split_longest(n3,P,lp,T,N,Bc,lb,Ph,Sh,L,Lx,af);
		else
			[Bc lb Sh] = split_boundary(-n3,Bc,lb,lp, Sh);
		end
	end
end % split_longest()

function [Bc lb Sh] = split_boundary(bdx,Bc,lb,pc,Sh)
	lb = lb+1;
	Bc(lb,1) = pc;
	Bc(lb,2) = Bc(bdx,2);
%	Bc(lb,3) = Bc(bdx,3);
	Bc(bdx,2) = pc;
	% push sides
	Sh.put(hashkey(Bc(bdx,1),Bc(bdx,2)),-bdx);
	Sh.put(hashkey(Bc(lb,1),Bc(lb,2)),-lb);
end

function [T Sh] = remove_edges(tdx,T,Sh)
	% remove old edge references
	% TODO, this has to go into split_edge and split_longest not into partition
	% avoid repitative removal of same entry
	% reason: old reference might be used, if one triangle is still ond and the other is new
	key = hashkey(T(tdx,1), T(tdx,2));
	val = Sh.remove(key)';
	if (length(val) > 2)
		if (tdx == val(1))
			Sh.put(key,val(3:4));
		else
			Sh.put(key,val(1:2));
		end
	else
		if (length(val) > 0 && tdx ~= val(1)) % already removed
			Sh.put(key,val(1:2));
		end
	end
	key = hashkey(T(tdx,1), T(tdx,3));
	val = Sh.remove(key)';
	if (length(val) > 2)
		if (tdx == val(1))
			Sh.put(key,val(3:4));
		else
			Sh.put(key,val(1:2));
		end
	else
		if (length(val) > 0 && tdx ~= val(1))
			Sh.put(key,val(1:2));
		end
	end
	key = hashkey(T(tdx,2), T(tdx,3));
	val = Sh.remove(key)';
	if (length(val) > 2)
		if (tdx == val(1))
			Sh.put(key,val(3:4));
		else
			Sh.put(key,val(1:2));
		end
	else
		if (length(val) > 0 && tdx ~= val(1))
			Sh.put(key,val(1:2));
		end
	end
end

function [T lt N Sh] = partition(tdx,T,lt,N,Ph,L,Lx,Sh)


	switch(Lx(tdx))
		case {1}
			p1 = T(tdx,1);
			p2 = T(tdx,2);
			p3 = T(tdx,3);
			l1=L(tdx,1);
			l2=L(tdx,2);
			l3=L(tdx,3);
		case {2}
			p1 = T(tdx,2);
			p2 = T(tdx,3);
			p3 = T(tdx,1);
			l1=L(tdx,2);
			l2=L(tdx,3);
			l3=L(tdx,1);
		case {3}
			p1 = T(tdx,3);
			p2 = T(tdx,1);
			p3 = T(tdx,2);
			l1=L(tdx,3);
			l2=L(tdx,2);
			l3=L(tdx,1);
	end

	% check, wether this triangle is obtuse
	cos_alpha = (l2^2 + l3^2 - l1^2)/(2*l2*l3);
	obtuse = (cos_alpha+sqrt(eps) < 0);

	% fetch new points
% TODO, use the T(4:6) and also preallocate
	%p12 = T(idx,6);
	%p13 = T(idx,5);
	%p23 = T(idx,4);
	p12 = Ph.get(hashkey(p1,p2));
	p13 = Ph.get(hashkey(p1,p3));
	p23 = Ph.get(hashkey(p2,p3));

	% determine, which sides in addition to 1-2 are split
	s = (~isempty(p23)) + 2*(~isempty(p13)) + 4*(~isempty(p12));
	% partition the triangle according to 4 different cases
	switch (s)
		case {0}
			% no side split
		case {1} % only longest side was split
			[T N Sh] = add_triangle(tdx,p1,p2,p23,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p1,p23,p3,T,N,Sh);
		case {3} % longest and right side were split
			[T N Sh] = add_triangle(tdx,p1,p2,p23,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p1,p23,p13,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p23,p3,p13,T,N,Sh);
		case {5} % longest and left side were split
			[T N Sh] = add_triangle(tdx,p1,p12,p23,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p2,p23,p12,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p1,p23,p3,T,N,Sh);
		case {7} % regular refinement, longest and both other sides were split
		if (~obtuse)
			[T N Sh] = add_triangle(tdx,p12,p13,p23,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p1,p12,p13,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p2,p12,p23,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p3,p13,p23,T,N,Sh);
		else
			[T N Sh] = add_triangle(tdx,p1,p12,p23,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p2,p23,p12,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p1,p23,p13,T,N,Sh);
			lt = lt+1;
			[T N Sh] = add_triangle(lt,p23,p3,p13,T,N,Sh);
		end
		otherwise
			error('here','longest side was not split');
	end
end % function partition

function [T N Sh] = add_triangle(tdx,p1,p2,p3,T,N,Sh)
	T(tdx,:) = [p1 p2 p3];
	% neighbour 1
	key = hashkey(p2,p3);
	val = Sh.remove(key);
	if (isempty(val))
		Sh.put(key,[tdx 1]);
	else
		if (val(1) > 0)
			% neighbour is a triangle
			N(tdx,1) = val(1);
			N(val(1),val(2)) = tdx;
		else
			% neighbour is a boundary
			N(tdx,1) = val(1);
		end
	end
	% neighbour 2
	key = hashkey(p1,p3);
	val = Sh.remove(key);
	if (isempty(val))
		Sh.put(key,[tdx 2]);
	else
		if (val(1) > 0)
			% neighbour is a triangle
			N(tdx,2) = val(1);
			N(val(1),val(2)) = tdx;
		else
			% neighbour is a boundary
			N(tdx,2) = val(1);
		end
	end
	% neighbour 3
	key = hashkey(p1,p2);
	val = Sh.remove(key);
	if (isempty(val))
		Sh.put(key,[tdx 3]);
	else
		if (val(1) > 0)
			% neighbour is a triangle
			N(tdx,3) = val(1);
			N(val(1),val(2)) = tdx;
		else
			% neighbour is a boundary
			N(tdx,3) = val(1);
		end
	end
end % add_triangle

