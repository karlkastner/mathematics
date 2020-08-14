% Sun  9 Aug 21:27:33 +08 2020
function y = mccormack_step(t,x,y,fun,dt)
	% dx same for both semi-steps
	dx = diff(x);
	yc = mid(y);

	f  = fun(y);
	fc = fun(yc);

	yc_ = yc - dt./dx*(f(2:end) - fc)
	y_  = inner2outer(yc_);
	f   = fun(y_)
	fc  = fun(yc_)
	yc  = 0.5*(yc + yc_) - 0.5*dt./dx*(fc - f(1:end-1))
	y   = inner2outer(yc)
end

