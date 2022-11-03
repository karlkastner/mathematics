% Tue  1 Nov 16:12:56 CET 2022
function w=tukeywin_man(x,p)
	w = zeros(size(x));
	% scale
	x = x - min(x);
	x = x/max(x);
	p = min(p,0.5);
	fdx = x<p & x > 0;
	w(fdx) = 0.5*(1.0 - cos(pi*x(fdx)/p));
	w(x>=p&x<=1-p) = 1;
	fdx = (x>1-p) & (x < 1);
	w(fdx) = 0.5*(1.0 - cos(pi*(x(fdx)-1)/p));
end

