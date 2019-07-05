% 2017-08-11 21:24:25.659767222 +0200
function s = sum_ij_partial_(rho,k,n,m)
	%N = 1:n-k;
	N = k+1:n;
	M = 1:m;
	%N = k+1:n;
	D = bsxfun(@minus,N',M);
	%D = bsxfun(@minus,N'+k,M);

%	N = k+1:n;
%	D = bsxfun(@minus,N',M);
	
	s = sum(sum(rho.^abs(D)));
end

