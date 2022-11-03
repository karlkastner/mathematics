% 2015-08-14 14:59:47.214583256 +0200
% Karl Kastner, Berlin

function test_sum_multivar
	rho1 = 0.25;
	rho2 = 0.5;
	s = [];
	n = 10;
	for idx=1:n
		s(idx,1) = sum_multivar(rho1,rho2,idx);
		s(idx,2) = sum_multivar_(rho1,rho2,idx);
	end
	s
end

function s = sum_multivar_(r1,r2,m)
	N = (0:m-1)';
	D = bsxfun(@minus,N,N');
	s = 0;
	O = ones(m);
	R = r1*triu(O,0) + r2*tril(O,0);
	s = sum(sum(R.^abs(D)));
end

