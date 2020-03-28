% Fri 17 Jan 12:33:40 +08 2020
function y = digamma_man(x)
	h = abs(x)*1e-5;
	d = (igamma(x+h,0)-igamma(x-h,0))/(2*h);
	y = d/igamma(x,0);
end
