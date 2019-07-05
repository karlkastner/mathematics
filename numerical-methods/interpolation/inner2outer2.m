% Fri  1 Jun 11:12:03 CEST 2018
%% interpolate from element (segment) centres to edge points
%
function x = inner2outer2(x)
	n  = size(x,1);
	x0 = x(1,:);
%	x0 = 0.5*(x(end,:)-x(1,:));
	L  = x(end,:)-x(1,:);
	x0 = x0-1/2*L/(n-1);
	x = bsxfun(@plus,x0,(0:n)'/(n-1)*L);
%	2*(n/2+1)/n*linspace(-1/2,1/2,n+1)'*L);
end

