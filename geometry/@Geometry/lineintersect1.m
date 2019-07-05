% 2015-10-28 19:40:15.756023461 +0100
%
%% intersect of two lines
function [flag s t p q den] = intersect(p1,p2,q1,q2)
	% determinant
	den = (p1(1)*q1(2) - p1(2)*q1(1) - p1(1)*q2(2) + p1(2)*q2(1) - p2(1)*q1(2) + p2(2)*q1(1) + p2(1)*q2(2) - p2(2)*q2(1));
	% normalised distance from intersection
	s  = (p1(1)*q1(2) - p1(2)*q1(1) - p1(1)*q2(2) + p1(2)*q2(1) + q1(1)*q2(2) - q1(2)*q2(1)) / den;
	t  = (-p1(1)*p2(2) + p1(2)*p2(1) + p1(1)*q1(2) - p1(2)*q1(1) - p2(1)*q1(2) + p2(2)*q1(1)) / den;
	% convexity (lines segments cut if convex)
	% this is strict, thus common end points are not regarded as intersection
	flag = (s > 0) & (s < 1) & (t > 0) & (t < 1);
	% intersection coordinates
	% p and q coincide only in case of intersection
	p(1) = (q1(1) + t*(q2(1)-q1(1)));
	p(2) = (q1(2) + t*(q2(2)-q1(2)));
	q(1) = (p1(1) + s*(p2(1)-p1(1)));
	q(2) = (p1(2) + s*(p2(2)-p1(2)));
end % intersect


