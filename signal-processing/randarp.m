% Fr 12. Feb 18:54:12 CET 2016
%
%% randomly generate the instance of an ar-p process
function [x D s] = randarp(sigma,r,n,m)
%	r = [r1 r2];
%	s2 =  (1-r(2))/(1+r(2))*1/((1-r(2))^2 - r(1)^2);
	%s = sqrt(s2);
	s = 1;
	R = randn(n,m);

	nr = length(r);
	D  = spdiags(ones(n,1)*[-r(end:-1:1), 1],-nr:0,n,n);
	x  = D \ R;
	% scale
	sc = std(flat(x));
	x  = sigma*bsxfun(@times,x,1./sc);
	s = 1./sc;

%	D = D/sqrt(s2);
%end
%	norm(x-X)
%	x = [x X];
end

