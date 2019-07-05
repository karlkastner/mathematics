% Di 1. Mär 14:04:36 CET 2016
%% simple quantile regression
function c = quatile_regression_simple(x,y,p)
	% iniial condition, ls line
	x = cvec(x);
	y = cvec(y);
	c0 = [x.^0 x] \ y;
%	c(1) = mean(y);
%	c(2) = 0;
	%tic()
	c = lsqnonlin(@(c) objective(x,y,c,p), c0);
	%toc()
	%tic()
	%c = fminsearch(@(c) objective(x,y,c,p), c0);
	%toc()
%	c(2) = lsqnonlin(@(c2) objective(x,y,[c(1);c2],p), c(2));
end

function f = objective(x,y,c,p)
	fdx = (y > c(1) + x*c(2));
	f = p*sum(abs(y(fdx) -(c(1) +  x(fdx)*c(2)))) + (1-p)*sum(abs(y(~fdx) - (c(1) + x(~fdx)*c(2))));
end

