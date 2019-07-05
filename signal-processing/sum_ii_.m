% 2017-08-11 16:31:42.923106775 +0200
function s = sum_ii_(rho,n)
	N = (1:n)';
	D = bsxfun(@minus,N,N');
	s = sum(sum(rho.^abs(D)));
end


