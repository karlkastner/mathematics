% 2015-06-23 12:27:24.983950647 +0200
% Karl Kastner, Berlin

function test_sum_ij

n = 10;
rho = 0.75;
s = [];
for m = 1:n
	s(m,1) = sum_ij(rho,n,m);
	s(m,2) = sum_ij_partial_(rho,m,n,n);
	%s(m,2) = sum_ij_(rho,n,m);
	%s(m,3) = doublesum_ij(rho,n,m);
	%s(m,3) = sum_ij_partial_(rho,m,n,n);
	%s(m,4) = (sum_ii(rho,n)-sum_ij_partial_(rho,n-m,n,n));
	%s(m,5) = (sum_ii(rho,n)-sum_ij_partial_(rho,n-m,n,n));
	%s(m,6) = (sum_ii(rho,n)-sum_ij(rho,n-m+1,n));
	%s(m,6) = sum_ij(rho,n,n-m);
	%s(m,5) = sum_ij_partial_(rho,n-m,n);
	%s(m,6) = sum_ij(rho,n-m,n);
end
'sum ij'
disp(s)
%s(:,3)-s(:,6)
'sum ii'	
[sum_ii(rho,n) sum_ii_(rho,n)]

'mu'
[sum_ii(rho,n)/n^2 mu2ar1(1,rho,n)]

n = 10
k = 6;
'sum i lag'
[sum_i_lag(rho,n,k) sum_i_lag_(rho,n,k)]

end


function s = sum_i_lag_(rho,n,k)
	N = (1:n)';
	D = N-k;
	s(1) = sum(rho.^abs(D));
end


