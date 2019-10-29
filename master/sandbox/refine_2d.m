% Tue May  1 19:15:48 MSK 2012
% Karl Kästner, Berlin
% 
% 
% 
% G: half split triangles information (0 == not half split, else index of other half)
% 
function [P T Bc N G] = refine_2d(P, T, Bc, N, M, G)
	import java.util.Hashtable;

	% get length
	lp = size(P,1);
	lt = size(T,1);
	lb = size(Bc,1);
	ld = 0;

	% allocate memory
	lp0 = 2*lp;
	lt0 = 4*lt;
	lb0 = 2*lb;
	P  = [P; zeros(lp,2)];
	T  = [T; zeros(3*lt,3)];
	Bc = [Bc; zeros(lb,3)];
	N  = [N; zeros(3*lt,3)];
	G  = [G; zeros(3*lt,1)];
	% effectively only 3/10 of the total number of triangles are at most merged
	D  = zeros(ceil(3*lt/10),1);

	% hanging nodes : hash of splitted sides, which are no boundary sides
	% key : lower_point_index, higher_point_index
	% value : new_point_index
	H = java.util.Hashtable(3*lt);

	% split marked elements into four parts
	for idx=1:length(M)
		[T lt P lp Bc lb N G H D ld] = split4(M(idx), T, lt, P, lp, Bc, lb, N, G, H, D, ld);
	end % for idx
	
	% remove hanging nodes
	while ( 1 ~= H.isEmpty() )
		% iterate over keys such that insertion and deletion is possible
		keySet = H.keySet();
		K = keySet.toArray();
		for kdx = 1:length(K)
			key = K(kdx);
			% remove hanging node from hash
			val = H.remove(key);
			% split into two, split2 might actualy split into four, to prevent degeneracy
			[T lt P lp Bc lb N G H D ld] = split2(val, T, lt, P, lp, Bc, lb, N, G, H, D, ld, key);
		end % for keys in hanging node list
	end  % while hanging nodes

	% remove merged triangle halves
	[T lt N G] = remove_merged(T, lt, N, G, D, ld);

	if (lp > lp0 || lt > lt0 || lb > lb0)
		warning('preallocation insufficient');
	end

	% reduce size to actual size
	P = P(1:lp,:);
	T = T(1:lt,:);
	N = N(1:lt,:);
	Bc = Bc(1:lb,:);
	G = G(1:lt,1);
end % fem_2d_refine

function [T lt P lp Bc lb N G H D ld] = split4(idx, T, lt, P, lp, Bc, lb, N, G, H, D, ld)
	% test if this is to be split but was already merged and 4-splitted
	% TODO dangerous ? should actually never happen
	if (0 == T(idx,1)) return; end

	% only split into four, if this is not a half triangle
	if (G(idx) == 0)
		% indices of the three points of the interior triangle
		T_ = zeros(1,3);
		% split each side
		for j1=1:3
			% jdx : index of opposit side (rotating index)
			% index-index of points bordering side
			j2 = mod(j1,3)+1;
			j3 = mod(j1+1,3)+1;
			% index of side points
			p2 = T(idx,j2);
			p3 = T(idx,j3);
			key = hashkey(p2,p3);
			% try to remove the entry from the hanging node list
			val = H.remove(key);
			% test if this side was already split for a neighbouring triangle
			if (isempty(val))
				% if not, create a new point (hanging node)
				lp = lp + 1;
				P(lp,:) = 0.5*(P(p2,:) + P(p3,:));
				T_(j1) = lp;
				% check if the new point is on the domain boundary
				if (N(idx,j1) < 0)
	 				% split the boundary segment in two parts
					Bc(-N(idx,j1),1:2) = [p2 lp];
					lb = lb + 1;
					Bc(lb,:) = [lp p3 Bc(-N(idx,j1),3)];
					% link from triangles to boundary segment
					 N(lt+j2,3) = N(idx,j1);
					 N(lt+j3,2) = -lb;
				else
					% not a boundary, neighbouring triangle exists
					% add hanging node
					if (0 == N(idx,j1))
						error('fem_2d_refine','here');
					end
					H.put(key, [lp N(idx,j1) lt+j2 lt+j3 p2]);
					% linking not yet possible, waiting until neighbour is split as well
				end
			else	% side was already split for neighbouring triangle
				pc = val(1);
				T_(j1) = pc;
				% check for 1-level recursion due to remerging
				% todo : not necessary, with small change in merge2
				val_ = H.get(hashkey(p2,pc));
				if (~isempty(val_))
					val_(2) = lt+j2;
					H.put(hashkey(p2,pc),val_);
					N(lt+j2, 3) = val_(5);
				else
					% link chidren of this triangle with its neighbours
					s2 = val(4);
					if (s2 > 0)
					if (p2 == T(s2,1))
						N(lt+j2, 3) = s2;
						N(s2,2) = lt+j2;
					elseif (p2 == T(s2,2))
						N(lt+j2, 3) = s2;
						N(s2,3) = lt+j2;
					elseif (p2 == T(G(s2,1)))
						N(lt+j2, 3) = G(s2);
						N(G(s2),2) = lt+j2;
					elseif (p2 == T(G(s2),2))
						N(lt+j2, 3) = G(s2);
						N(G(s2),3) = lt+j2;
					else
						error('fem_2d_refine','here');
					end
					end
				end % recursion not isempty val_
				% check for 1-level recursion due to remerging
				% todo : not necessary, with small change in merge2
				val_ = H.get(hashkey(p3,pc));
				if (~isempty(val_))
					val_(2) = lt+j3;
					H.put(hashkey(p3,pc),val_);
					N(lt+j3, 2) = val_(5);
				else
					% link chidren of this triangle with its neighbours
					% four cases
					s1 = val(3);
					if (s1 > 0)
					if (p3 == T(s1,1))
						N(lt+j3, 2) = s1;
						N(s1,3) = lt+j3;
					elseif (p3 == T(s1,3))
						N(lt+j3, 2) = s1;
						N(s1,2) = lt+j3;
					elseif (p3 == T(G(s1),1))
						N(lt+j3, 2) = G(s1);
						N(G(s1),3) = lt+j3;
					elseif (p3 == T(G(s1),3))
						N(lt+j3, 2) = G(s1);
						N(G(s1),2) = lt+j3;
					else
						error('fem_2d_refine','here');
					end
					end
				end % recursion not isempty val_
			end % not isempty val
		end % for j1
		% old points (exterior)
		a = T(idx,1); 	b = T(idx,2);	c = T(idx,3);
	
		% new points (interior)
		A = T_(1); B = T_(2); C = T_(3);
	
		% add exterior triangles
		T(lt+1,:) = [a C B];
		T(lt+2,:) = [b A C];
		T(lt+3,:) = [c B A];
	
		% overwrite interior triangle
		T(idx,:) =  [A B C];
	
		% link exterior triangles with interior triangle
		N(lt+1,1) = idx;
		N(lt+2,1) = idx;
		N(lt+3,1) = idx;
		
		% link interior triangle with exterior triangles
		N(idx,:) = lt+[1 2 3];
		lt = lt+3;
	else % this was a half triangle, merge and split it into four
		[T lt P lp Bc lb N G H D ld] = merge2(idx, T, lt, P, lp, Bc, lb, N, G, H, D, ld);
		[T lt P lp Bc lb N G H D ld] = split4(idx, T, lt, P, lp, Bc, lb, N, G, H, D, ld);
	end
end % function split4

% triangle halving to remove hanging nodes
% quarters triangle in case of degeneracy
function [T lt P lp Bc lb N G H D ld] = split2(val, T, lt, P, lp, Bc, lb, N, G, H, D, ld, key)
	% allready implicitely splitted
	if (isempty(val))
		return;
	end
	idx = val(2);

	if (0 == idx)
		errplot(T,P,B);
	 	error('fem_2d_refine','here');
	end

	% check that this triangle was not yet split into halves
	if (0 == G(idx))
		pc  = val(1);
		t1  = val(3);
		t2  = val(4);
		pf  = val(5);
		% find the corner point opposit the hanging node
		for k1=1:3
			k3 = mod(k1+1,3)+1;
			p3 = T(idx,k3);
			if (pf == p3)
				lt = lt+1;
				k2 = mod(k1,3)+1;
				p1 = T(idx,k1);
				p2 = T(idx,k2);
				n2 = N(idx,k2);
				n3 = N(idx,k3);

				% check if partner was halved, four cases
				if (p3 == T(t1,1))
					N(t1,3) = lt;
				elseif (p3 == T(t1,3))
					N(t1,2) = lt;
				elseif (G(t1) == 0)
					error('','');
				elseif (p3 == T(G(t1),1))
					t1 = G(t1);
					N(t1,3) = lt;
				elseif (p3 == T(G(t1),3))
					t1 = G(t1);
					N(t1,2) = lt;
				else
					error('fem_2d_refine','here');
				end
				if (p2 == T(t2,1))
					N(t2,2) = idx;
				elseif (p2 == T(t2,2))
					N(t2,3) = idx;
				elseif (G(t2) == 0)
					errplot(T,lt,P,lp,Bc,lb,N,G);
					error('fem_2d_refine','here');
				elseif (p2 == T(G(t2),1))
					t2 = G(t2);
					N(t2,2) = idx;
				elseif (p2 == T(G(t2),2))
					t2 = G(t2);
					N(t2,3) = idx;
				else
					error('fem_2d_refine','here');
				end
				
				% order similar to 4-split corner triangles
				T(lt,:)  = [ p3, p1, pc];
				N(lt,:)  = [idx, t1, n2];
				T(idx,:) = [ p2, pc, p1];
				N(idx,:) = [ lt, n3, t2];

				G(idx) = lt;
				G(lt) = idx;
				% if neighbour is not a boundary
				if (n2 > 0)
					key = hashkey( p1, p3 );
					val = H.get(key);
					if (~isempty(val))
						val(2) = lt;
						H.put(key,val);
					else % empty
						for l1=1:3
							if (N(n2,l1) == idx)
								N(n2,l1) = lt;
								break;
							end % if ldx
							if (3 == l1) error('fem_2d_refine','here'); end
						end % for
					end % not isempty
				end % if n3 > 0
				break;
			end % if jdx
			if (3 == k1)
					disp(idx)
					subplot(2,2,3); pdx = find(T(1:lt,1) > 0); disp([ pdx T(pdx,:) N(pdx,:)])
					display_2d(P(1:lp,:),T,Bc(1:lb,:),7,ones(size(pdx)),pdx); axis equal; axis tight; set(gca,'xtick',[]); set(gca,'ytick',[]);
					axis equal; axis tight; title('after remove'); drawnow; pause(2);
					error('fem_2d_refine','here');
			end
		end % for k1
	else % triangle was alread split into halves, merge and split into four
		% repush this side to be splitted in halves
		H.put(key,val);
		% merge halves and split into four
		[T lt P lp Bc lb N G H D ld parent] = merge2(idx, T, lt, P, lp, Bc, lb, N, G, H, D, ld);
		[T lt P lp Bc lb N G H D ld] = split4(parent, T, lt, P, lp, Bc, lb, N, G, H, D, ld);
	subplot(2,2,3)
	tdx=find(T(:,1));
	bdx=find(Bc(:,1));
	mesh = Mesh(P, T(tdx,:), Bc(bdx,:));
	display_2d(mesh, 7, [], [], 'EdgeColor', 'k');
%                K = H.keySet().toArray();                for idx=1:length(K)                          disp( [key (H.get(K(idx)))']);                 end

%	idx
%	display_2d(mesh, 7, [], [], 'EdgeColor', 'k');
%	pause
	end % merge-split
end % function split2

% merges two halv triangles and split the merged triangle into four parts
% idx becomes the merged triangle
function [T lt P lp Bc lb N G H D ld parent] = merge2(idx, T, lt, P, lp, Bc, lb, N, G, H, D, ld)
	% triangle index of other half
	gdx = G(idx);

	% find point in triangle gdx opposit triangle idx
	for k1=1:3
		if (idx == N(gdx,k1))
			break;
		end % if idx == N(gdx,k1)
		if (3 == k1) errplot(T,lt,P,lp,Bc,lb,N,G);
				disp([idx G(idx) T(idx,:) N(idx,:); gdx G(gdx) T(gdx,:) N(gdx,:)]);
				 error('fem_2d_refine','here - do not flip!');
		end
	end % for k1
	% find point in triangle idx opposit triangle gdx
	for l1=1:3
		if (gdx == N(idx,l1))
			break;
		end % if gdx==N(idx,l1)
		if (3 == l1) error('fem_2d_refine','here'); end
	end % for l1
	% mark triangles not to be halves any more
	G(idx) = 0;
	G(gdx) = 0;

	% extract indices
	q1 = T(gdx,k1);
	p1 = T(idx,l1);
	l2 = mod(l1,3)+1;
	l3 = mod(l1+1,3)+1;
	p2 = T(idx,l2);
	p3 = T(idx,l3);

	% two scenarios
	A = [1 P(p1,:); 1 P(q1,:); 1 P(p2,:)];
	B = [1 P(p1,:); 1 P(p2,:); 1 P(p3,:)];
	if (abs(det(A)/det(B)) > 1e-12)
		tmp = idx; idx=gdx; gdx=tmp;
		tmp = k1;    k1=l1; l1=tmp;
	end

	% properties of first triangle
	k2 = mod(k1,3)+1;
	k3 = mod(k1+1,3)+1;
	m1 = N(gdx,k1);
	m2 = N(gdx,k2);
	m3 = N(gdx,k3);
	q1 = T(gdx,k1);
	q2 = T(gdx,k2);
	q3 = T(gdx,k3);

	% properties of second triangle
	l2 = mod(l1,3)+1;
	l3 = mod(l1+1,3)+1;
	p1 = T(idx,l1);
	p2 = T(idx,l2);
	p3 = T(idx,l3);
	n3 = N(idx,l3);

	T(idx,l2) = q1;
	N(idx,l1) = m3;

	% check if m3 is not a boundary
	if (m3 > 0)
		% check if side is to be splitted
		key = hashkey(q1,q2);
		val = H.get(key);
		if (~isempty(val))
			val(2) = idx;
			H.put(key,val);
		else
		% unsplitted neighbour exists
			for mdx=1:3
				if (N(m3,mdx) == gdx)
					N(m3,mdx) = idx;
					break;
				end % if N(m1,mdx) == gdx
				if (3 == mdx) error('fem_2d_refine','here'); end
			end % for mdx
		end % if isempty
	end % if m1 > 0
	% check if m2 is not a boundary
	if (m2 > 0)
		% check if this side is to be splitted
		key = hashkey(q1,q3);
		val = H.get(key);
		if (~isempty(val))
			val(2) = idx;
			H.put(key,val);
			m2_ = 0;
		else
		% unsplitted neighbour exists
			m2_ = m2;
			for mdx=1:3
				if (N(m2,mdx) == gdx)
					N(m2,mdx) = idx;
					break;
				end % if N(m1,mdx) == gdx
				if (3 == mdx) error('fem_2d_refine','here'); end
			end % for mdx
		end % if isempty
	end % if m2 > 0
	% check if side n3 is to be splitted
	val = H.get(hashkey(p1,q2));
	if (~isempty(val))
		n3_ = 0;
	else
		n3_ = n3;
	end

	% add a hanging node
	val = [p2 idx m2_ n3_ q1];
	H.put(hashkey(p1,q1),val);
	% reset the dead triangle
	T(gdx,:) = 0;
	N(gdx,:) = 0;
	% mark dead triangle to be deleted
	ld = ld+1;
	D(ld) = gdx;
	parent = idx;
end % merge2

function [T lt N G] = remove_merged(T, lt, N, G, D, ld)
	% sort to avoid cases like remove last, remove last but two -> last but one lost
	D = sort(D,1,'descend');
	for idx=1:ld
		ddx = D(idx);
		if (ddx < lt)
		T(ddx,:) = T(lt,:);
		N(ddx,:) = N(lt,:);
		if (G(lt) > 0)
			G(G(lt)) = ddx;
			G(ddx) = G(lt);
		end
		for kdx=1:3
			% not for boundaries
			if (N(ddx,kdx) > 0)
				for jdx=1:3
					if (N(N(ddx,kdx),jdx) == lt)
						N(N(ddx,kdx),jdx) = ddx;
						break;
					end
					if (3 == jdx) error('fem_2d_refine','here'); end
				end % for jdx
			end % ifN(ddx,kdx) > 0
		end % for kdx
		end
		lt = lt-1;
	end % for idx
end % remove_merged

function errplot(T,lt,P,lp,Bc,lb,N,G)
		%subplot(2,2,3);
		pdx = find(T(1:lt,1) > 0); disp('after remove'), disp([ pdx T(pdx,:) N(pdx,:) G(pdx)]);
		%display_2d(P(1:lp,:),T,Bc(1:lb,:),7,ones(size(pdx)),pdx); axis equal; axis tight; set(gca,'xtick',[]); set(gca,'ytick',[]);
		mesh = Mesh(P(1:lp,:),T(1:lt,:),Bc(1:lb,:));
		clf
		display_2d(mesh,0,ones(size(pdx)),pdx); axis equal; axis tight; set(gca,'xtick',[]); set(gca,'ytick',[]);
		axis equal; axis tight; title('after remove'); drawnow; pause(2);
end

