% Sun Jun 15 18:28:16 WIB 2014
% Karl Kastner, Berlin

%% transform points from cartesian into streamwise coordinates
%%
%% NOTE : prefer the java version, this has some problems with round off
% TODO this seems to suffer from numerical instability,
%      so maybe it is better to translate A into 0 (as now)
%	and then rotate AB to be parallel to the x axis
%	and than simply find delta S and N as X and Y coordinates of C
% TODO : how to better deal with multiple sections and bifurcations ?
% todo 10s progress indicator
function [S N nn] = xy2sn(cX,cY,cS,X,Y)
	linear = 1;
	% for each point in X, find closest centreline point and distance to it
	S = zeros(size(X));
	N = zeros(size(X));
	% find nearest neighbours on the centreline for each input point
	nn = knnsearch([cX cY],[X Y],'K',1);

	nok = 0;
	for idx=1:length(X)
		% distance of current sampled point to all centreline points
		%d = (cX - X(idx)).^2 + (cY - Y(idx)).^2;
		% find index of closesest centreline point
		%[mv mx1] = min(d);
		mx1 = nn(idx);
	
		if (linear)
			% determine the second neighbour
			% TODO this should not be second closest but
			%      the neighbour that makes the transformation function convex,
			%      e.g. [X Y] = p [cX1 + (1-p) cX2, 0 <= p <= 1
			% TODO, this may crash at the end !!!
			%d2left  = dist2([cX(mx1) cY(mx1)],[cX(mx1-1) cY(mx1-1)]);
			%d2right = dist2([cX(mx1) cY(mx1)],[cX(mx1+1) cY(mx1+1)]);
			%d2left  = dist2([X(idx) Y(idx)],[cX(mx1-1) cY(mx1-1)]);
			%d2right = dist2([X(idx) Y(idx)],[cX(mx1+1) cY(mx1+1)]);

			if (mx1-1 < 1)
				d2left = inf;
			else
				d2left = dist2([X(idx) Y(idx)],[cX(mx1-1) cY(mx1-1)]);
			end
			if (mx1+1 > length(cX))
				d2right = inf;
			else
				d2right = dist2([X(idx) Y(idx)],[cX(mx1+1) cY(mx1+1)]);
			end

			if (d2left < d2right)
				mx2 = mx1-1;
				mx0 = mx1+1;
				mx3 = mx2-1;
			else
				mx2 = mx1+1;
				
				mx0 = mx1-1;
				mx3 = mx2+1;
			end
			% corner vertices of the triangle
			A = [cX(mx1); cY(mx1)];
			B = [cX(mx2); cY(mx2)];
			C = [X(idx); Y(idx)];
			% make A the local origin
			B = B-A;
			C = C-A;
			A = [0;0];

			% squared side length of triangle
			AB = B-A;
			AC = C-A;
			BC = C-B;
			aa = BC'*BC;
			bb = AC'*AC;
			cc = AB'*AB;
			% height of triangle
			n = sqrt(2*(aa*bb + bb*cc + aa*cc) - (aa*aa + bb*bb + cc*cc))/(2*sqrt(cc));

			% distance to foot point from mx1
			% ds = bb/sqrt(cc) applies only true to right angled triangles
			% and misses direction, so min ||C-P||, s.t. P = A + alpha AB yields:
			%ds = (C'*AB - A'*AB)/sqrt(cc);
			%P = A+ds*AB/sqrt(cc);
			% s.t. P = pA + (1-p)B
			p = (B'*B - A'*B + A'*C - B'*C)/(A'*A - 2*A'*B + B'*B );
			P = p*A + (1-p)*B;
			PC = C-P;
			s = p*cS(mx1) + (1-p)*cS(mx2);
			% n is unsigned, so determine, whether it is left or right of AB
			n = sign(PC(2)*AB(1)*(cS(mx2)-cS(mx1)))*n;
		
			% also not obvious from the formula, the sign is random
	%		P1 = A+ds*AB/sqrt(cc);
	%		P2 = A-ds*AB/sqrt(cc);
			%ds = sign(AC(1)*AB(1)*(cS(mx2)-cS(mx1)))*ds;
	%		[ (C-P1)'*(C-P1) (C-P2)'*(C-P2)]
	%		ds_ = ds;
	%		if ( (C-P1)'*(C-P1) < (C-P2)'*(C-P2) )
	%			ds = abs(ds);
	%		else
	%			ds = -abs(ds);
	%		end
	%		[ds_ ds]
	%		ds = -sign(cS(mx2)-cS(mx1))*ds;
	%		s = cS(mx1)+ds;

			%if (cS(mx1) < cS(mx2))
			%	s = cS(mx1)+ds;
			%else
			%	s = cS(mx1)-ds;
			%end
		else
			% constant interpolation (warning: leads to severe aliasing)
			s = cS(mx1);
			% TODO sign
			n = sqrt(mv);
		end % if not linear
		S(idx) = s;
		N(idx) = n;
	end % for idx
end % xy2sn

