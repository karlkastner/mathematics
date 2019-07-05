% Do 11. Feb 18:51:53 CET 2016
% Karl Kastner, Berlin
%
%% random number of laplace distribution
%
function r = laplacernd(n,m)
	r1 = exprnd(1,n,m);
	r2 = exprnd(1,n,m);
	r = r1-r2;
end

