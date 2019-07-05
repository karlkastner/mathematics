% 2016-03-03 19:22:07.947910317 +0100
% Karl Kastner, Berlin
%
%% pdf of the bivariate normal distribution
function z = normpdf2(XY,mu,S)
	dxy = bsxfun(@minus,XY,mu);
	Si  = inv(S);
	dS  = det(S);
	z   = 1./sqrt(2*pi*dS)*exp(-0.5*sum(dxy.*(Si*dxy')',2));
end
