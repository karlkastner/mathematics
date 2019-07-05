% Thu Dec 11 15:09:08 CET 2014
% Karl Kastner, Berlin
%
%% fit linear function ||a b x = y||_L2 by the method of moments
%% y+eps = alpha + beta*x
%
% TODO this is not really robust, median should be chosen for alpha
function [c, s2res] = mreg(x,y)
	n = length(x);
	mux = mean(x);
	muy = mean(y);
	% variances
	s2x = var(x);
	s2y = var(y);
	% covariance
	C = cov(x,y);
	cxy = C(1,2);
	% estimate
	beta = cxy/s2x;
	alpha = muy - beta*mux;
	% residual variance
	s2res = s2y - beta*beta*s2x;
		
	% normalise
	s2res = s2res*(n-2)/(n-1);
	
	% solution vector for [x.^0 x.^1] c = y + eps
	c = [alpha; beta];
end

