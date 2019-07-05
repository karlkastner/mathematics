% 2013-06-23 11:53:40.000000000 +0200
%
%% maximum value and i-j index for matrix
function [v, x, y] = max2d(X)
	[v,ind]=max(X(:));
	[y,x]  = ind2sub(size(X),ind);
end

