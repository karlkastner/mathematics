% 2019-06-27 19:01:29.521493388 +0200
function s = sum_ij_(rho,n,m)
	N = (1:n)';
	M = (1:m)';
	D = bsxfun(@minus,N,M');
	s = sum(sum(rho.^abs(D)));
end

