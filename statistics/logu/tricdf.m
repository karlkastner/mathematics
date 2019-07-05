% Thu 22 Mar 17:27:07 CET 2018
% Karl Kastner, Berlin
%% cumulative distribution of the log-triangular distribution
function F = logtricdf(a,b,c,x)
	F      = (x-a).^2/((c-a)*(b-a));
	F(x<a) = 0;
	fdx    = x>b;
	F(fdx) = 1-(c-x(fdx)).^2/((c-a)*(c-b));	
	F(x>c) = 1;
end

