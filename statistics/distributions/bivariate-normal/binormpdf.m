% Tue Nov 18 23:11:21 CET 2014
% Karl Kastner, Berlin
%
% probability density of a bi-modal normal distribution
%
%pdf = @(x,p,mu1,mu2,s1,s2) p*normpdf(x,mu1,s1) + (1-p)*normpdf(x,mu2,s2);
function p = binormpdf(x,p,mu1,mu2,sd1,sd2)
	p  =      par(1)*normpdf(x,par(2),par(4)) ...
	     + (1-par(1))*normpdf(x,par(3),par(5));
end

