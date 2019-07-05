% Di 22. Dez 15:01:05 CET 2015
% Karl Kastner, Berlin
%
%% intersections of two polygons
function [X0 Y0 FDX DEN] = polyxpoly(X,Y)
	X   = rvec(X);
	Y   = rvec(Y);
	X0  = [];
	Y0  = [];
	n   = length(X);
	FDX = [];
	DEN = [];
	for idx=1:n-1
		p1 = [X(idx);Y(idx)];
		p2 = [X(idx+1);Y(idx+1)];
		q1 = [X(idx+1:n-1);
                      Y(idx+1:n-1)];
		q2 = [X(idx+2:n);
                      Y(idx+2:n)];
		[f s t p q den] = Geometry.lineintersect(p1,p2,q1,q2);
		fdx = cvec(find(f));
		if (length(fdx) > 0)
			FDX = [FDX; ...
                               [idx*ones(size(fdx)), idx+fdx]];
			X0 = [X0; p(1,fdx)'];
			Y0 = [Y0; p(2,fdx)'];
			DEN = [DEN; den(fdx)];

		end
	end
end % polyxpoly

