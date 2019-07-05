% Fri Jul  4 19:01:06 WIB 2014
% Karl Kastner, Berlin
%
%% convert cartesian to streamwise coordiantes
%
function [S N] = xy2sn2(cX,cY,cS,X,Y)

	% for each point in X, find closest centreline point and distance to it
	S = zeros(size(X));
	N = zeros(size(X));
	% find nearest neighbours on the centreline for each input point
	nn = knnsearch([cX cY],[X Y],'K',1);

	for idx=1:length(X)
		% get the index of closest centreline point
		mx1 = nn(idx);
	
		% determine second neighbour
		% TODO this should not be second closest but
		%      the neighbour that makes the transformation function convex,
		%      e.g. [X Y] = p [cX1 + (1-p) cX2, 0 <= p <= 1
		if (1 == mx1)
			mx2 = mx1+1;
		elseif (mx1 == length(cX))
			mx2 = mx1-1;
		else
			d2left  = Geometry.distance2([X(idx) Y(idx)],[cX(mx1-1) cY(mx1-1)]);
			d2right = Geometry.distance2([X(idx) Y(idx)],[cX(mx1+1) cY(mx1+1)]);
			if (d2left < d2right)
				mx2 = mx1-1;
			else
				mx2 = mx1+1;
			end
		end
		% corner vertices of the triangle
		A = [cX(mx1); cY(mx1)];
		B = [cX(mx2); cY(mx2)];
		C = [X(idx); Y(idx)];

		% translace coordinates such that A is the origin
		B = B-A; 
		C = C-A;
		% A := 0
		% rotate coordinates, such that AB = (B-0) becomes oriented along the x-axis
		% C := R^-1 C
		hyp = sqrt(B'*B);
		sina = B(2)/hyp;
		cosa = B(1)/hyp;
		Ri = [ cosa sina
	              -sina cosa];
		C = Ri*C;
		% now Cx is the dS and Cy is N
		if (cS(mx1) < cS(mx2))
			S(idx) = cS(mx1) + C(1);
			N(idx) = C(2);
		else
			S(idx) = cS(mx1) - C(1);
			N(idx) = -C(2);
		end
%		S_(idx,:) = [cS(mx1) S(idx) cS(mx2)]-cS(mx1);
%		catch e
%			e
%			S(idx) = NaN;
%			N(idx) = NaN;
%		end
	end % for idx
%	close all
%	plot([S_]);
%	pause

end % xy2sn

