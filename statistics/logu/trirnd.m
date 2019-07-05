% Thu 22 Mar 17:14:50 CET 2018
% Karl Kastner, Berlin
%
%% random numbers of the triangular distribution
function y = trirnd(a,b,c,m,n)
	x = rand(m,n);
	y = triinv(a,b,c,x);
end

