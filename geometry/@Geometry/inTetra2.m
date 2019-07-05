% Thu 11 Jan 19:09:28 CET 2018
% Karl Kastner, Berlin
%
%% flag points contained in tetrahedron
function [flag, c] = inTetra2(P1,P2,P3,P4,P0)

	c = Geometry.tobarycentric3(P1, P2, P3, P4, P0);

	tol  = eps^(1/3);
	% tol = eps^0.25;
	flag = prod(c >= -tol & c <= 1+tol);

	flag = logical(flag);
	flag = cvec(flag);
end

