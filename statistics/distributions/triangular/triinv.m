% Thu 22 Mar 17:15:40 CET 2018
%
%% inverse of the triangular distribution
function [x] = triinv(a,b,c,F)
	F0     = (b-a)/(c-a);
%	F0     = (c-b)/(c-a);
	x      = a + sqrt(F*(c-a)*(b-a));
	fdx    = F>F0;
	x(fdx) = c - sqrt((1-F(fdx))*(c-a)*(c-b));
	x(F<0) = NaN;
	x(F>1) = NaN;
	% TODO, switch case
	% x = exp((a + sqrt(F)))*(a-c)*(a-b);
end

