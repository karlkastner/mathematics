% 2014-06-29 15:34:22.826392717 +0200
% Karl Kastner, Berlin
%
%% intersect of two lines
% TODO normalisation, why is there just one denominator ?
% what is s, what t?
% function [flag, s, t, p, q, den] = intersect(p1,p2,q1,q2)
% flag : true if intersecting
% s,t  : relative distance on p and q to intersection point
% p, q : closest point to other section
function [flag, s, t, p, q, den] = intersect(p1,p2,q1,q2)
	np = size(p1,2);
	nq = size(q1,2);
	P1X = repmat(p1(1,:),nq,1);
	P1Y = repmat(p1(2,:),nq,1);
	P2X = repmat(p2(1,:),nq,1);
	P2Y = repmat(p2(2,:),nq,1);
	Q1X = repmat(q1(1,:)',1,np);
	Q1Y = repmat(q1(2,:)',1,np);
	Q2X = repmat(q2(1,:)',1,np);
	Q2Y = repmat(q2(2,:)',1,np);

	% determinant
if (0)
	den =   P1X.*Q1Y - P1Y.*Q1X - P1X.*Q2Y + P1Y.*Q2X ...
              - P2X.*Q1Y + P2Y.*Q1X + P2X.*Q2Y - P2Y.*Q2X;

	% normalised distance from intersection
	s  = (  P1X.*Q1Y - P1Y.*Q1X - P1X.*Q2Y + P1Y.*Q2X + Q1X.*Q2Y - Q1Y.*Q2X ) ./ den;
	t  = (- P1X.*P2Y + P1Y.*P2X + P1X.*Q1Y - P1Y.*Q1X - P2X.*Q1Y + P2Y.*Q1X ) ./ den;
else
	% backup
	Q1X_ = Q1X;
	Q1Y_ = Q1Y;

	% shift to avoid round off error
	P2X = P2X - P1X;
	Q1X = Q1X - P1X;
	Q2X = Q2X - P1X;
	% P1X = 0.*P1X;
	P2Y = P2Y - P1Y;
	Q1Y = Q1Y - P1Y;
	Q2Y = Q2Y - P1Y;
	% P1Y = 0.*P1Y;

	% because P1X and P1Y are zero:
	den =  - P2X.*Q1Y + P2Y.*Q1X + P2X.*Q2Y - P2Y.*Q2X;

	% normalised distance from intersection
	s  = (    Q1X.*Q2Y - Q1Y.*Q2X ) ./ den;
	t  = (  -(P2X.*Q1Y) + P2Y.*Q1X ) ./ den;
end

	% convexity (lines segments cut if convex)
	% this is strict, thus common end points are not regarded as intersection
	delta = sqrt(eps);
	flag = (s > delta) & (s < 1-delta) & (t > delta) & (t < 1-delta);
	% intersection coordinates
	% p and q coincide only in case of intersection
	if (nargout() > 3)
		p(1,:,:) = bsxfun(@plus, Q1X_, bsxfun(@times, t, (q2(1,:)-q1(1,:))'));
		p(2,:,:) = bsxfun(@plus, Q1Y_, bsxfun(@times, t, (q2(2,:)-q1(2,:))'));
	if (nargout() > 4)
		q(1,:,:) = bsxfun(@plus, P1X, bsxfun(@times, s, (p2(1,:)-p1(1,:))));
		q(2,:,:) = bsxfun(@plus, P1Y, bsxfun(@times, s, (p2(2,:)-p1(2,:))));
	end
	end
end % intersect

