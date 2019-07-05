% Thu 22 Mar 09:54:53 CET 2018
%
%% random numbers following a log-uniform distribution
%
% function x = logurnd(a,b,n,m)
function y = logurnd(a,b,m,n)
	x = rand(m,n);
	%y = a*exp(log(b/a)*x);
	%y = a*(b/a).^x;
	y = loguinv(a,b,x);
end

