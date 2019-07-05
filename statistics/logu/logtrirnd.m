% Fri 23 Mar 13:36:26 CET 2018
%
% random numers of the logarithmic triangular distribution
function x = logtrirnd(a,b,c,m,n)
	F = rand(m,n);
	x = logtriinv(a,b,c,F);
end

