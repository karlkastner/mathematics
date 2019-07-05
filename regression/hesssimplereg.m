% 2016-02-25 13:57:00.687966775 +0100
%
%% hessian, gradient and objective function value of the simple regression
%% rhs = p(1) + p(2) x + eps
%
% TODO get residual vectos
function [H, g, f] = hesssimplereg(x,rhs,p)
	n      = length(x);
%	H      = zeros(2,2);
	sx     = sum(x);
	H      = 2*[ n, sx;
		    sx, sum(x.^2)];
%	g      = zeros(2,1);
	res = p(1) + p(2)*x - rhs;
	g  = 2*[sum(res);
	        sum(x.*(res)];
	f  = sum(res.^2);
end

