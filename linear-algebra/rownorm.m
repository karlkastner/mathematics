% 2017-11-14 16:03:41.896408949 +0100
% norms of rows of x
function x = rownorm(x)
	x = sqrt(sum(x.^2,2));
end
