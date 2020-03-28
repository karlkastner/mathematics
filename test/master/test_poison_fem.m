N=2.^(0:8);

for idx=1:length(N)
	n = N(idx);
	[A B] = poisson_fem(n);
	L = chol(B);
	A=inv(L')*A*inv(L);
	err(idx,1) = eigs(A,[],1,'SM')+pi^2;
	A = A + 1/12*A^2/(n+1)^2;
	err(idx,2) = eigs(A,[],1,'SM')+pi^2;
	 %AA = A + 1/12*A^2/(n+1); BB = B + 12*B^2*(n+1); eigs(AA,BB,1,'SM')
end
loglog(N,abs(err))

