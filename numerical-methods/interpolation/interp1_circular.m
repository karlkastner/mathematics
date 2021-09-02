function yi = interp1_circular(x,y,xi,varargin)
	x = [2*x(1)-x(2);
             cvec(x);
	     2*x(end)-x(end-1)];
	y = [2*y(1)-y(2);
             cvec(y);
	     2*y(end)-y(end-1)];
	yi = interp1(x,y,xi,varargin{:});
end
