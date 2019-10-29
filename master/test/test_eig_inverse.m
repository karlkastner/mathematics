% Sun Dec  2 15:14:03 MSK 2012
% Karl Kästner, Berlin

function test_eig_inverse()
	n = 100;
	A = rand(n);
	A=A*A';
	
	[v e] = eig_inverse(A);
	e = eig(A);
	e=min(diag(e));

	semilogy([R abs(L-e)])
end % test_eig_inverse()

