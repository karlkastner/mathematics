% Thu 22 Jun 10:43:25 CEST 2017
% determine if c is located on the 1d line segment
function [flag c] = onLine(X,x0)
	c = (x0 - X(:,1))./(X(:,2)-X(:,1));
	flag = (c >= 0) & (c <= 1);
end
