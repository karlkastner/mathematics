% Tue Nov 18 23:12:21 CET 2014
% Karl Kastner,Berlin
%
%% bio-modal gaussian distribution
% pdf = @(x,p,mu1,mu2,s1,s2) p*normpdf(x,mu1,s1) + (1-p)*normpdf(x,mu2,s2);
%
function p = binormcdf(x,par)
	p  =       par(1)*normcdf(x,par(2),par(4)) ...
             + (1-par(1))*normcdf(x,par(3),par(5));
end

