% 2015-11-02 18:56:47.018003787 +0100
% Karl Kastner, Berlin
%
%% flag points contained in triangle
function [flag c] = inTriangle(P1,P2,P3,P0)

	c  = Geometry.tobarycentric2(P1,P2,P3,P0);
	%c = Geometry.tobarycentric2(X,Y,X0,Y0);

	tol  = eps^(1/3);
	% tol = eps^0.25;
	flag = prod(c >= -tol & c <= 1+tol);

	flag = logical(flag);
	flag = cvec(flag);
end

