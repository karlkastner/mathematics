% triangle halving to remove hanging nodes
% quarters triangle in case of degeneracy
function [T lt P lp Bc lb N Hh] = split2_final(val, T, lt, P, lp, Bc, lb, N, Hh)
		if (isempty(val))
			return
		end
		pc  = val(1);
		idx = val(2);
		t1  = val(3);
		t2  = val(4);
		jdx = val(5);

		% find the corner point opposit the hanging node
		for k1=1:3
			if (jdx == N(idx,k1))
				k2 = mod(k1,3)+1;
				k3 = mod(k1+1,3)+1;
				p1 = T(idx,k1); p2 = T(idx,k2); p3=T(idx,k3);
				n2 = N(idx,k2); n3 = N(idx,k3);
				A = [ 	1 P(T(idx,1),:);
					1 P(T(idx,2),:);
					1 P(T(idx,3),:) ];
					% angle on opposit side is to small (30°)
				if (abs(det(A))/(norm(P(p2,:) - P(p1,:))*norm(P(p3,:) - P(p1,:))) < 0.5) % && ...
					 %norm(P(p2,:)-P(p1,:)) < norm(P(p1,:) - P(p3,:)) ) 
					if ( norm(P(p2,:)-P(p1,:)) > norm(P(p1,:) - P(p3,:)) ) 
 						% split side p2-p1
 						% check, if hanging node exists
 						val = Hh.remove(hashkey(p1,p2));
 						if (isempty(val))
 							lp = lp+1;
 							pd = lp;
 							P(lp,:) = 0.5*(P(p2,:)+P(p1,:));
 							% check for boundary
 							if ( n3 < 0)
 								lb = lb+1;
 								Bc(-n3,1:2) = [p2 pd];
 								Bc(lb,:) = [pd p1 Bc(-n3,3)];
 								s1 = n3;
 								s2 = -lb;
 							else
 								% add hanging node
 								% TODO the 5th value can be avoided
 								val = [pd n3 idx lt+1 idx]; % lt+2
 								Hh.put(hashkey(p1,p2),val);
 								s1 = 0; %inf;
 								s2 = 0; %-inf;
 							end
 						else
 							pd = val(1);
 							s1 = val(3);
 							s2 = val(4);
 							N(s1,3) = lt+1;
 							N(s2,2) = idx;
 						end % isempty
 						T(idx,:)   = [p1 pd p3];
 						T(lt+1,:)  = [p2 pc pd];
 						T(lt+2,:)  = [p3 pd pc];
 						N( idx,:)  = [lt+2 n2  s2];
 						N(lt+1,:)  = [lt+2 s1  t2];
 						N(lt+2,:)  = [lt+1 t1 idx];
 						N(t1,3)    = lt+2;
 						N(t2,2)    = lt+1;
 						lt = lt+2;
				 % }
				else
						% split side p3-p1
						% check, if hanging node exists
						val = Hh.remove(hashkey(p1,p3));
						if (isempty(val))
							lp = lp+1;
							pd = lp;
							P(lp,:) = 0.5*(P(p3,:)+P(p1,:));
							% check for boundary
							if ( n2 < 0 )
								lb = lb+1;
								Bc(-n2,1:2) = [p3 pd];
								Bc(lb,:) = [pd p1 Bc(-n2,3)];
								b  =  n2;
								a  = -lb;
							else
								% add hanging node
								val = [pd n2 lt+1 idx idx];
								Hh.put(hashkey(p1,p3),val);
								a  = 0;
								b  = 0;
							end
						else
							pd = val(1);
							a  = val(3);
							b  = val(4);
							N( a,3) = idx; %lt+1;
							N( b,2) = lt+1; %idx;
						end % isempty
						T(idx,:)   = [p1   p2  pd];
						N( idx,:)  = [lt+2  a  n3];
						T(lt+1,:)  = [p3   pd  pc];
						N(lt+1,:)  = [lt+2 t1   b];
						T(lt+2,:)  = [p2   pc  pd];
						N(lt+2,:)  = [lt+1 idx t2];
						N(t1,3)    = lt+1;
						N(t2,2)    = lt+2;
						'Molch'
						%[p1 p2 p3 pd]
						[idx T(idx,:) N(idx,:)]
						[lt+1 T(lt+1,:) N(lt+1,:)]
						[lt+2 T(lt+2,:) N(lt+2,:)]
						lt = lt+2;
					end % split p1-p3
% }
				else % angle sufficiently large				
				lt = lt+1;
% %{
% 				lt = lt+1;
% 				% do not rotate !
% 				 %T(idx,:) = [ p1,  pc, p3];
% 				 %N(idx,:) = [ t1,  n2, lt];
% 				 T(idx,k2) = pc;
% 				 N(idx,k1) = t1;
% 				 N(idx,k3) = lt;
% 				N(t1,3) = idx;
% 				T(lt,:) = [ p1,  p2, pc];
% 				N(lt,:) = [ t2, idx, n3];
% 				N(t2,2) = lt;
% 				'honk 3'
% 				[idx T(idx,:)]
% 				[lt T(lt,:)]
% 				% if neighbour is not a boundary
% 				if (n3 > 0)
% 					key = hashkey( p1, p2 );
% 					val = Hh.get(key);
% 					if (~isempty(val))
% 						val(2) = lt;
% 						Hh.put(key,val);
% 					else % empty
% 						for l1=1:3
% 							if (N(n3,l1) == idx)
% 								N(n3,l1) = lt;
% 								break;
% 							end % if ldx
% 						end % for
% 					end % not isempty
% 				end % if n3 > 0
% %}
			if (0 == n2)
				% do not rotate !
				%T(idx,:) = [ p1,  pc, p3];
				%N(idx,:) = [ t1,  n2, lt];
				 T(idx,k2) = pc;
				 N(idx,k1) = t1;
				 N(idx,k3) = lt;
				N(t1,3) = idx;
				T(lt,:) = [ p1,  p2, pc];
				N(lt,:) = [ t2, idx, n3];
				N(t2,2) = lt;
				'honk 2'
				[idx T(idx,:) N(idx,:)]
				[lt T(lt,:) N(lt,:)]
				% if neighbour is not a boundary
				if (n3 > 0)
					key = hashkey( p1, p2 );
					val = Hh.get(key);
					if (~isempty(val))
						val(2) = lt;
						Hh.put(key,val);
					else % empty
						for l1=1:3
							if (N(n3,l1) == idx)
								N(n3,l1) = lt;
								break;
							end % if ldx
						end % for
					end % not isempty
				end % if n3 > 0
 			else % 0 == n3
				T(idx,k3) = pc;
				N(idx,k1) = t2;
				N(idx,k2) = lt;
				T(lt,:)   = [p1 pc  p3];
				N(lt,:)   = [t1 n2 idx];
				N(t1,3)   = lt;
				N(t2,2)   = idx;
				'honk 3'
				[idx T(idx,:) N(idx,:)]
				[lt T(lt,:) N(lt,:)]
				if (n2 > 0)
					key = hashkey(p1, p3);
					val = Hh.get(key);
					if (~isempty(val))
						val(2) = lt;
						Hh.put(key,val);
						'val'
					else
						for l1=1:3
							if (N(n2,l1) == idx)
								N(n2,l1) = lt;
								l1
								break;
							end
						end % for l1
					end % if isempty
				end % if n2 > 0
			end % if 0 == n3
				end % angle sufficiently large
				break;
			end % if jdx
		end % for k1
end % function split2



% triangle halving to remove hanging nodes
% quarters triangle in case of degeneracy
function [T lt P lp Bc lb N Hh] = split2(val, T, lt, P, lp, Bc, lb, N, Hh)
		pc  = val(1);
		idx = val(2);
		t1  = val(3);
		t2  = val(4);
		jdx = val(5);
		% find the corner point opposit the hanging node
		for k1=1:3
			if (jdx == N(idx,k1))
				% check if this triangle and a neighbour are already a half split
				k2 = mod(k1,3)+1;
				k3 = mod(k1+1,3)+1;
				n2 = N(idx,k2);
				n3 = N(idx,k3);
				p1 = T(idx,k1); p2 = T(idx,k2); p3=T(idx,k3);
				if (n2 > 0)
					for ldx=1:3
						if (T(n2,ldx) == p1)
							break;
						end
					end
					n21 = N(n2,ldx);
					n22 = N(n2,mod(ldx,3)+1);
					p23 = T(n2,mod(ldx+1,3)+1);
					A2 = [1 P(p1,:);
 					       1 P(p3,:);
					       1 P(p23,:)];
				else
					A2 = inf;
				end % n2 > 0
				if (n3 > 0)
					for ldx=1:3
						if (T(n3,ldx) == p1)
							break;
						end
					end
					p22 = T(n3,mod(ldx,3)+1);
					A3 = [ 1 P(p2,:);
 					       1 P(p3,:);
					       1 P(p22,:)];
				else
					A3 = inf;
				end % n3 > 0
				if (abs(det(A2)) < 1e-12)
					% save indix
					sdx = n2;
					% remerge
					T(idx,3) = p24;
					N(idx,3) = n22;
					if (n22 > 0)
						for ldx=1:3
							if (n2 == N(n22, ldx))
								N(n22,ldx) = idx;
								break;
							end
						end % for ldx
					end % for n22
					% add hanging node created by merging
					Hh.put(hashkey(p2,p4), p3, n21, t2, idx);
					% split into 4
					[T lt P lp Bc lb N Hh] = split4(idx, T, lt, P, lp, Bc, lb, N, Hh);
					% resolve the hanging node
					T(n2,:)    = [ pc p3    T(idx,k3)];
					N(n2,:)    = [idx lt+k2 t1];
					N(idx,k2)  = n2;
					N(t1,3)    = n2;
					T(lt+k2,2) = pc;
					N(lt+k2,1) = n2;
					N(lt+k2,3) = t2;
					N(t2,2) = lt+k2;
				 %elseif (abs(det(A24)) < 1e-12)
				 %	TODO
				else
				lt = lt+1;
				T(idx,:) = [ p1,  pc, p3];
				N(idx,:) = [ t1,  n2, lt];
				N(t1,3) = idx;
				T(lt,:) = [ p1,  p2, pc];
				N(lt,:) = [ t2, idx, n3];
				N(t2,2) = lt;
				% if neighbour is not a boundary
				if (n3 > 0)
					key = hashkey( p1, p2 );
					val = Hh.get(key);
					if (~isempty(val))
						val(2) = lt;
						Hh.put(key,val);
					else % empty
						for l1=1:3
							if (N(n3,l1) == idx)
								N(n3,l1) = lt;
							break;
							end % if ldx
						end % for
					end % not isempty
				end % if n3 > 0
				end % if two-split
				break;
			end % if jdx
		end % for k1
end % function split2

%{
						% check if right also must be splitted
						val = Hh.remove(hashkey(p1,p3));
						if (~isempty(val))
							pe = val(1);
							r1 = val(3);
							r2 = val(4);
							T(lt+1,:) = [pd p3 pe];
							T(idx,:)  = [p1 pd pe];
							N(lt+1,:) = [r2 idx lt];
							N(lt,3)   = lt+1;
							N(idx,:)  = [lt+1 r1 s2];
							N(r1,3)   = idx;
							N(r2,2)   = lt+1;
							lt = lt+1;
						end
%}
%{

						T(lt+1,:)  = [p2 pc pd];
						N(lt+1,:)  = [lt+2 s1  t2];
						T(lt+2,:)  = [p3 pd pc];
						N(t1,3)    = lt+2;
						N(t2,2)    = lt+1;
						% check if right also must be splitted
						val = Hh.remove(hashkey(p1,p3));
						if (~isempty(val))
							pe = val(1);
							r1 = val(3);
							r2 = val(4);
							N(lt+2,:) = [lt+1 t1 lt+3];
							T(idx,:)  = [p1 pd pe];
							N(idx,:)  = [lt+3 r1 s2];
							T(lt+3,:) = [pd p3 pe];
							N(lt+3,:) = [r2 idx lt+2];
							N(r1,3)   = idx;
							N(r2,2)   = lt+3;
							lt = lt+3;
						else
							N(lt+2,:) = [lt+1 t1 idx];
							T( idx,:)  = [p1   pd  p3];
							N( idx,:)  = [lt+2 n2  s2];
							lt = lt+2; 
						end
%}

	%{
					else	% split side p1-p3
						% check, if hanging node exists
						val = Hh.remove(hashkey(p1,p3));
						if (isempty(val))
							lp = lp+1;
							pd = lp;
							P(lp,:) = 0.5*(P(p3,:)+P(p1,:));
							% add hanging node
									% the 5th value can be avoided
							Hh.put(hashkey(p1,p3),[pd N(idx,k3) lt+1 idx lt+2]);
							s1 = 0;
							s2 = 0;
						else
							pd = val(1);
							s1 = val(3);
							s2 = val(4);
							N(s1,2) = idx;
							N(s2,3) = lt+1;
						end % isempty
						T(idx,:)   = [p1   p2 pd];
						T(lt+1,:)  = [pd   pc p3];
						T(lt+2,:)  = [pd   p2 pc];
						N(idx,1:2) = [s1 lt+2];
						N(lt+1,:)  = [lt+2 s2  t2];
						N(lt+2,:)  = [lt+1 idx t1];
						N(t2,3)    = lt+2;
						N(t1,2)    = lt+1;
						lt = lt+2;
					end % split p1-p3
% }
					else
%}
% triangle halving to remove hanging nodes
% quarters triangle in case of degeneracy
function [T lt P lp Bc lb N Hh] = split2_4(val, T, lt, P, lp, Bc, lb, N, Hh)
		if (isempty(val))
			return
		end
		pc  = val(1);
		idx = val(2);
		t1  = val(3);
		t2  = val(4);
		jdx = val(5);

		% find the corner point opposit the hanging node
		for k1=1:3
			if (jdx == N(idx,k1))
				k2 = mod(k1,3)+1;
				k3 = mod(k1+1,3)+1;
				p1 = T(idx,k1); p2 = T(idx,k2); p3=T(idx,k3);
				n2 = N(idx,k2); n3 = N(idx,k3);
				A = [ 	1 P(T(idx,1),:);
					1 P(T(idx,2),:);
					1 P(T(idx,3),:) ];
					% angle on opposit side is to small (30°)
				if (abs(det(A))/(norm(P(p2,:) - P(p1,:))*norm(P(p3,:) - P(p1,:))) < 1/sqrt(2))
					Hh.put(hashkey(p2,p3),val);
					[T lt P lp Bc lb N Hh] = split4(idx, T, lt, P, lp, Bc, lb, N, Hh);
				else % angle sufficiently large				
				lt = lt+1;
			if (0 == n2)
				% do not rotate !
				%T(idx,:) = [ p1,  pc, p3];
				%N(idx,:) = [ t1,  n2, lt];
				 T(idx,k2) = pc;
				 N(idx,k1) = t1;
				 N(idx,k3) = lt;
				N(t1,3) = idx;
				T(lt,:) = [ p1,  p2, pc];
				N(lt,:) = [ t2, idx, n3];
				N(t2,2) = lt;
				% if neighbour is not a boundary
				if (n3 > 0)
					key = hashkey( p1, p2 );
					val = Hh.get(key);
					if (~isempty(val))
						val(2) = lt;
						Hh.put(key,val);
					else % empty
						for l1=1:3
							if (N(n3,l1) == idx)
								N(n3,l1) = lt;
								break;
							end % if ldx
						end % for
					end % not isempty
				end % if n3 > 0
 			else % 0 == n3
				T(idx,k3) = pc;
				N(idx,k1) = t2;
				N(idx,k2) = lt;
				T(lt,:)   = [p1 pc  p3];
				N(lt,:)   = [t1 n2 idx];
				N(t1,3)   = lt;
				N(t2,2)   = idx;
				if (n2 > 0)
					key = hashkey(p1, p3);
					val = Hh.get(key);
					if (~isempty(val))
						val(2) = lt;
						Hh.put(key,val);
					else
						for l1=1:3
							if (N(n2,l1) == idx)
								N(n2,l1) = lt;
								break;
							end
						end % for l1
					end % if isempty
				end % if n2 > 0
			end % if 0 == n3
			end % if angle sufficiently large
			break;
			end % if jdx
		end % for k1
end % function split2

