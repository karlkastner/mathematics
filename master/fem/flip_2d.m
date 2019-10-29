% 2012 May  6 17:54 (MSK)
% Karl Kästner, Berlin

function [P T Bc N] = flip(M, P, T, Bc, N, G)
	for mdx=1:length(M)
		idx = M(mdx);
		% only try to flip if this triangle is not green!
		if (0 ~= G(idx)) continue; end;
		for j1=1:3
			n1 = N(idx,j1);
			% onle flip, if the neihbour is not a boundary
			if (n1 > 0)
			% only try to flip if this triangle is not green!
			if (0 ~= G(n1)) continue; end
			j2 = mod(j1,3) + 1;
			j3 = mod(j1+1,3) + 1;
			n2 = N(idx,j2);
			%n3 = N(idx,j3);
			p1 = T(idx,j1);
			p2 = T(idx,j2);
			p3 = T(idx,j3);
%			[idx n1 n2]
			% find opposit point
			for k1=1:3
				if (N(n1,k1) == idx)
					p4 = T(n1,k1);
					k2 = mod(k1,3) + 1;
					k3 = mod(k1+1,3) + 1;
					m2 = N(n1,k2);
					break;
				end
			end
			% check if opposit < common
			if ( norm(P(p1,:) - P(p4,:)) < norm(P(T(idx,j2),:) - P(T(idx,j3),:))-1e-12 )
				% neither p2 nor p3 have to be convex combinations
				% check that the fourth point is a non-convex combination
				A = [ 1 P(p1,:);
				      1 P(p2,:);
				      1 P(p4,:) ]';
				if (abs(det(A)) > 1e-12)
					c3 = A \ [1 P(p3,:)]';
				else
					c3 = 0;
				end
				A = [ 1 P(p1,:);
				      1 P(p3,:);
				      1 P(p4,:) ]';
				if (abs(det(A)) > 1e-12)
					c2 = A \ [1 P(p2,:)]';
				else
					c2 = 0;
				end
				% sum(c) == 1, always, not necessary to check
				if ( (sum(c2 > 0)) < 3 && sum(c3 > 0) < 3 &&  sum( abs(c2) < 1e-12) == 0 && sum(abs(c3) < 1e-12) == 0 )
					% not a convex combination, flip possible
					T(idx,j3) = p4;
					N(idx,j1) = m2; 
					N(idx,j2) = n1;
					T(n1,k3)  = p1;
					N(n1,k2)  = idx;  
					N(n1,k1)  = n2;
					if (n2 > 0)
					for ldx=1:3
						if (N(n2,ldx) == idx)
							N(n2,ldx) = n1;
							break;
						end
					end
					end
					if (m2 > 0)
					for ldx=1:3
						if (N(m2,ldx) == n1)
							N(m2,ldx) = idx;
							break;
						end
					end
					end
				end
			end
			end
		end % for j1
	end % for mdx
end % function flip

