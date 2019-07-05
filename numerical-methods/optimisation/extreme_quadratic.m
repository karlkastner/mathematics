% 2017-03-30 10:58:00.490966245 +0200
function [x0 y0] = extreme_quadratic(x,y)
	% shift
	xc = x(:,2);
	yc = y(:,2);
	x = bsxfun(@minus,x,xc);
	y = bsxfun(@minus,y,yc);
	c = ((y(:,1)).*x(:,3)    - y(:,3).*x(:,1))./(x(:,1).^2.*x(:,3) - x(:,3).^2.*x(:,1));
	b = ((y(:,1)).*x(:,3).^2 - y(:,3).*x(:,1).^2)./(x(:,1).*x(:,3).^2 - x(:,3).*x(:,1).^2);
	x0 = -0.5*b./c;
	y0 = b.*x0 + c.*x0.^2;
	
	% shift back
	x0 = x0 + xc;
	y0 = y0 + yc;
end


