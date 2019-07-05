% Thu 22 Mar 10:53:02 CET 2018
%
%% probability density of the logarithmic uniform distribution
function F = logucdf(a,b,x)
	F = log(x/a)/log(b/a);
	F(x<a) = 0;
	F(x>b) = 1;
end

