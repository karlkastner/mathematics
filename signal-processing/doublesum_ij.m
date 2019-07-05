% 2015-08-05 17:01:21.980006350 +0200
% Karl Kastner, Berlin
%
%% double sum of r^i
function s = doublesum_ij(rho,n,k)
	s(1) = (n-k)*(1+rho)/(1-rho) ...
		- 1/(1-rho)^2*(rho - rho^(n-k+1) + rho^(k+1) - rho^(n+1));
	N = (1:n)';
	M = (1:n-k);
	D = bsxfun(@minus,N,M);
%	s(2) = sum(sum(rho.^(abs(D))));
%	s(3) = sum(sum(rho.^(abs(D-k))));
%	s(3) = 0;
%	l = 0;
%	r = 0;
%	for idx=1:n-k
%		%s(2) = s(2)+sum_ij(rho,idx,n);
%		s(3) = s(3) - 1/(1-rho)*(rho^(idx)+rho^(n-idx+1));
%		l = l - 1/(1-rho)*rho^idx;
%		r = r - 1/(1-rho)*rho^(n-idx+1);
%	end
%	[l -1/(1-rho)^2*(rho-rho^(n-k+1))]
%	[r]
%	s(3) = s(3)+(n-k)*(1+rho)/(1-rho);
%	s(4) = (n-k)*(1+rho)/(1-rho)+r+l;
end

