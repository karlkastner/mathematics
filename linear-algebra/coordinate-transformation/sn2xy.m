% Mo 8. Sep 18:34:59 CEST 2014
% Karl Kastner, Berlin
%
%% convert sn to xy coordinates
function [X, Y, sobj] = sn2xy(cX,cY,cS,S,N,sobj)
	k = length(S);
	X = zeros(k,1);
	Y = zeros(k,1);
	% create search object
	if (nargin() < 6 || isempty(sobj))
		sobj = createns(cS);
	end
	M1 = knnsearch(sobj,S);
	for idx=1:k
		m1 = M1(idx);
		% choose neighbour so that S is a convex combination
		if (1 == m1)
			% left boundary
			m2 = m1+1;
		elseif (length(cS) == m1)
			% right boundary
			m2 = m1-1;
		elseif((cS(m1) - S(idx))*(cS(m1+1) - S(idx)) < 0)
			m2 = m1+1;
		else
			m2 = m1-1;
		end
		ds = (S(idx) - cS(m1));
		dn = N(idx); % N of centre is zero
		if (cS(m2) < cS(m1))
			ds = -ds;
			dn = -dn;
		end
		% direction of centre line
		dx = cX(m2) - cX(m1);
		dy = cY(m2) - cY(m1);
		% TODO the hypothenuse equals abs(S2-S1) if S was not calculated along an arc
		hyp = hypot(dx,dy); 
		% perform rotation and translation
		c = dx/hyp;
		s = dy/hyp;
		X(idx) = c*ds - s*dn + cX(m1);
		Y(idx) = s*ds + c*dn + cY(m1);
	end % for idx
end % function sn2xy

