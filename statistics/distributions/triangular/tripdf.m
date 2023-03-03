% Thu 22 Mar 17:20:40 CET 2018
% Karl Kastner, Berlin
%
%% probability density of the triangular distribution
function f = tripdf(a,b,c,x)
	f      = 2*(x-a)/((c-a)*(b-a));
	f(x<a) = 0;
	fdx = x>b;
	f(fdx) = 2*(c-x(fdx))/((c-a)*(c-b));
	f(x>c) = 0;
end

