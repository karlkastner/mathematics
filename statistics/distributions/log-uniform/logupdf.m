% Thu 22 Mar 10:53:42 CET 2018
% Karl Kastner, Berlin
%
%% pdf of the log uniform distribution
function f = logupdf(a,b,x)
	f = 1./(x*log(b/a));
	f(x<a) = 0;
	f(x>b) = 0;
end

