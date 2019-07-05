% Mon 26 Sep 13:34:06 CEST 2016
% Karl Kastner, Berlin
%
%% edge lenght of triangle
function [d2 d] = tri_side_length(X,Y)
	dX = (left(X)-right(X));
	dY = (left(Y)-right(Y));
	d2 = dX.^2 + dY.^2;
	if (nargout() > 1)
		d = sqrt(d2);
	end
end

