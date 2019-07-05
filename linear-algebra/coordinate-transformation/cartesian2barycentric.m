% Fri 23 Sep 20:14:34 CEST 2016
%
%% cartesian to barycentric coordinates
function [pq] = cartesian2barycentric(xy0,xy)
	A = [xy0(1,1:2)-xy0(1,3);
	     xy0(2,1:2)-xy0(2,3)];
	b = bsxfun(@minus,xy,xy0(:,3));

	pq = A \ b;
end

