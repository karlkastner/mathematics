% 2016-05-09 17:48:36.721621179 +0200
%% slope by linear regression
function [slope] = slope(x,y,W)
	if (isvector(x))
		x = cvec(x);
		y = cvec(y);
	end
	if (nargin() < 3)
		W     = [];
	end

	x     = bsxfun(@minus,x,mean(x));
	order = 1;
	param = PolyOLS.fit_(x,y,W,order);
	slope = param(2);
end

