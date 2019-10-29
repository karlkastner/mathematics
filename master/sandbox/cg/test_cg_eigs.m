% n >= 9, no convergence due to round off
% largest eigenvalus converge at first

	n=5;
	A = rand(n);
	A=A*A';
	errmat(E, flipud(eig(A)))

	m = 3; n=5;
	A = - laplacian_2d(m,n);
	n = m*n;

	b=rand(n,1);
	[x nerrs a b_] = cg(A,b)
	y=A \ b
	norm(x-y)/norm(x)
	p = cg_coef_to_poly(a,b_)
	R=zeros(n);
	for m=1:n;
		p = cg_coef_to_poly(a(1:m),b_(1:m));
		R(1:m,m) = roots(p);
	end
	[R flipud(sort(eig(full(A))))]
	errmat(R,eig(full(A)))

