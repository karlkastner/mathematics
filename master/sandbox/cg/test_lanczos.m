	E=[];
	m = 4; n=5;
	A = - laplacian_2d(m,n);
	[a b T] = lanczos(A);
for idx=1:n*m
	E(1:idx,idx) = flipud(sort(eig(full(T(1:idx,1:idx)))));
end
	[E flipud(eig(A))]
M=errmat(E, flipud(eig(A)))
plot(M(:,end-5:end-1))

