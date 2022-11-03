% 2021-06-21 15:47:09.925953141 +0200
function [win,r] = window2d(f,r)
	n = length(f);
	% TODO, fourier axis
	x = 2*(0:n-1)/(n-1)-1;
	r = hypot(x,x');
	win = 0.5*(r<=1).*(1+cos(pi*r));
	win = win/mean(win(:));
end

