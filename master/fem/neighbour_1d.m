% Tue Jun 19 15:48:22 MSK 2012
% Karl Kästner, Berlin

function N = neighbour_1d(P,T)
	lt = size(T,1);
	lp = size(P,1);
	% find element neighbours
	N = zeros(lt,2);
	Q = zeros(lp,1);
	for idx=1:lt
		for jdx=1:2
			p = T(idx,jdx);
			if (0 == Q(p))
				Q(p) = idx;
			else
				if (0 == N(idx,1)) N(idx,1) = Q(p); else N(idx,1) = Q(p); end
				if (0 == N(Q(p),1)) N(Q(p),1) = idx; else N(Q(p),2) = idx; end
			end
		end % for jdx
	end % for idx
end % neighbour_1d

