% Sat 10 Sep 12:50:00 CEST 2016
%
%% orthocentre of triangle
function [x0 y0] = orthocentre(x,y)
	xl=left(x);
	xr=right(x);
	yl=left(y);
	yr=right(y);
	% orthogonal direction
	d = [yr-yl;
             xl-xr];
%	plot([x-d(1,:);x+d(1,:)],[y-d(2,:);y+d(2,:)])
	z = [0;0];
	b = -[x(1)-x(2);
	      y(1)-y(2);
      	      x(1)-x(3);
      	      y(1)-y(3);
	      x(2)-x(3);
	      y(2)-y(3)];
	A = [d(:,1),-d(:,2),z;
             d(:,1),z,-d(:,3);
             z,d(:,2),-d(:,3)];
	c = A \ b;
	% for numerical stability, we evaluate p at all points and average
	x0 = x + d(1,:).*c.';
	y0 = y + d(2,:).*c.';
	x0 = mean(x0);
	y0 = mean(y0);
end
