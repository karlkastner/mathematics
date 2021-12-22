% 2021-09-24 22:48:34.731934312 +0200 high_pass_simple.m
function y = high_pass_simple(x,a)
	% y(idx) = a*y(idx-1) + a*(x(idx)-x(idx-1));
	y(:,2) = filter(a*[1,-1],[1,-a],x);
end
