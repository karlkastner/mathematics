% 2016-03-30 15:01:02.377775205 +0200
% norms of rows of x
function [x nrow] = normrows(x)
	nrow = sqrt(sum(x.^2,2));
	x = bsxfun(@times,x,1./nrow);
end

