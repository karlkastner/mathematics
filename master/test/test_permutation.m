N=(1:30)';
M=[]

for idx=1:length(N);
	n = N(idx)
	I=speye(n);
	A=spdiags(ones(n,1)*[1 -2 1], -1:1, n,n);
	A=kron(kron(A,I),I)+kron(kron(I,A),I)+kron(kron(I,I),A);
	M(idx,1) = nnz(A);
	tic()
	M(idx,2) = nnz(chol(-A));
	T(idx,1) = toc();
	tic()
	p = symamd(A);
	M(idx,3) = nnz(chol(-A(p,p)))
	T(idx,2) = toc();
end
%S = [N.^4 N.^3 N.^2 N N.^0];
S = [ N.^5 N.^4 N.^3 N.^2 N.^1];
S \ M
subplot(1,2,1)
loglog(N,M)
subplot(1,2,2)
loglog(N,T)
M(:,2)./M(:,3)

