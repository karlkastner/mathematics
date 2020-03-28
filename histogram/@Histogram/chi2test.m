% 2014-11-22 21:50:50.098573016 +0100
% Karl Kastner, Berlin
function [alpha obj] = chi2test(obj,N1,N2)
%	h = h1.*h2;
%	n = (n1*n2)/(n1+n2);
%	chi2 = sum( (n-n).^2/njk
%	N1 = p1*n1;
%	N2
%	N1
%	N2
	n1   = sum(N1);
	n2   = sum(N2);
	n    = n1+n2;
	N1s  = (N1+N2)*n1/n;
	N2s  = (N1+N2)*n2/n;
	fdx = (N1+N2) > 0;
	chi2 = ((N1-N1s).^2)./N1s + ((N2-N2s).^2)./N2s
	chi2  = sum(chi2(fdx))
	dof  = sum(fdx)-1;
%	dof  = length(N1)-1;
	alpha = 1-chi2cdf(chi2,dof);
	%alpha = chi2inv(chi2,dof%);

%	N1s = round(N1s)
%	N2s = round(N2s)
%	(N1-N1s).^2
end

