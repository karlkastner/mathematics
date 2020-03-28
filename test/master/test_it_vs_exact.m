%Tue Jun  5 10:59:13 MSK 2012
%Karl Kästner, Berlin
function test_it_vs_chol(pflag)
	path(path,'../jacobi-davidson');
	path(path,'../lanczos');
	path(path,'../');
	path(path,'/home/pia/Desktop/cosse/kth/dn2230-numalg/hw1/');
%	subplot(2,2,1)
%	test_solvers();
	figure(2); clf;
	subplot(2,2,1)
	test_eigs(2,1);
	subplot(2,2,2)
	test_eigs(2,10);
	subplot(2,2,3)
	test_eigs(3,1);
	subplot(2,2,4)
	test_eigs(3,10);

	if (nargin()>0 && pflag)
		print -depsc eigensolver-poisson-run-times-i.eps
	end

function x = afunc_minres(A, L, x, tol, maxit)
	[x flag] = minres(A, x, tol, maxit, L',L);
end

function x = afunc_chol(L, Lt, p, x)
	x(p) = L \(Lt\x(p));
end

function x = afunc_lu(L, U, P, Q, D, x)
        x = Q*(U \ (L \ (P*(D\x))));
end

function test_solvers()
	tol = sqrt(eps);
	n0 = 30;

	% minres vs chol 1D
	n = n0^3;
	A = -spdiags(ones(n,1)*[1 -2 1], -1:1, n, n);
	b = -rand(n,1);
	opts.maxit = n;
	tic();
	L = chol(A);
	x = L\(L'\b);
	t(1) = toc()
	tic();
	% stagnates
	%x_ = minres(A,b,tol,opts.maxit);
	x_=x;
	t(2) = toc()
	
	% minres vs chol 2D
	n = round(n0^(3/2));
	A = -spdiags(ones(n,1)*[1 -2 1], -1:1, n, n);
	I = speye(n);
	A = kron(A,I) + kron(I,A);
	b = -rand(n^2,1);
	opts.maxit = n^2;

	tic();
	[L U pp qq R] = lu(A);
	tic()
        x = qq*(U \ (L \ (pp*(R\b))));
	toc()
	t(2)=toc()
	tic();
	p = symamd(A);
	L = chol(A(p,p));
	I = speye(size(A));
	P = I(p,p);
	Q(p,p) = I;
	Lt = L';
	tic()
	x(p) = L\(Lt\b(p));
	%x(p) = L\(L'\b(p));
	%x(p) = Q*(L\(L'\(P*b)));
	toc()
	t(3) = toc()
	% no p
	tic()
	[x_ flag] = minres(A,b,tol,opts.maxit,[],[]);
	t(4) = toc()
	tic();
	L = ichol(A);
	[x_ flag] = minres(A,b,tol,opts.maxit,L',L);
	toc()
	t(5) = toc()

	% minres vs chol 3D
	n = n0;
	A = -spdiags(ones(n,1)*[1 -2 1], -1:1, n, n);
	I = speye(n);
	A = kron(kron(A,I),I) + kron(kron(I,A),I) + kron(kron(I,I),A);
	b = -rand(n^3,1);
	opts.maxit = n^3;
	tic();
	p = symamd(A);
	L = chol(A(p,p));
	x(p) = L\(L'\b(p));
	t(6) = toc()
	toc()
	tic();
	L = ichol(A);
	[x_ flag] = minres(A,b,tol,opts.maxit,[],[]);
	t(7) = toc()
	tic();
	L = ichol(A);
	[x_ flag] = minres(A,b,tol,opts.maxit,L,L');
	t(8) = toc()
	norm(x-x_)
	bar(t);
	title(['Iterative vs Exact Solver ' num2str(n0) '^3 grid points']);
	set(gca,'xticklabel',{'Chol 1D', 'Minres 1D', 'Chol 2D', 'Minres 2D', 'PMinres 2D', 'Chol 3D', 'Minres 3D', 'PMinres 3D'});
end

function test_eigs(d, k)
	tol = sqrt(eps);
	opts.issym = 1;
	%opts.isreal = 0;
	opts.tol=tol;
	opts.disp=0;
	% minres vs chol 3D
	n = round(15^(3/d));
	opts.maxit = n^d;
	%opts.v0 = ones(n^d,1) + 0.1*rand(n^d,1);
	%opts.v0 = rand(n^d,1);
	A = -(n+1)^2*spdiags(ones(n,1)*[1 -2 1], -1:1, n, n);
	I = speye(n);
	switch(d)
		case {2}
			A = kron(A,I) + kron(I,A);
		case {3}
			A = kron(kron(A,I),I) + kron(kron(I,A),I) + kron(kron(I,I),A);
	end

%	tic();
%	jdopts.maxit = n^d;
%	eigs_lanczos(A, k, n^d, 'SM')
%	t(1) = toc()
	tic();
	jdopts.maxit = n^d;
	jdqr(A, k, 'SM', jdopts)
	t(1) = toc()
	
	tic();
	eigs(A,[], k,'SM',opts)
	t(2) = toc()
	% lu
	tic();
	[L U P Q D] = lu(A);
	eigs(@(x) afunc_lu(L,U,P,Q,D,x), n^d, 1, 'SM')
	t(3) = toc()
	tic();
	p = symamd(A);
	L = chol(A(p,p));
	Lt = L';
	%I = speye(n^d); P = I(p,p); Q(p,p) = I;
	eigs(@(x) afunc_chol(L,Lt,p,x), n^d, 1, 'SM')
	t(4) = toc()
	tic();
	L = ichol(A);
	eigs(@(x) afunc_minres(A,L,x,tol,n^d), n^d, k, 'SM', opts)
	t(5) = toc()
	bar(t);
	title(['Eigensolver Runtime for Poisson Equation in ' num2str(d) 'D with ' num2str(n) '^' num2str(d) ' Points and ' num2str(k) ' computed eigenvalues']);
	%set(gca,'xticklabel',{'Chol 1D', 'Minres 1D', 'Chol 2D', 'Minres 2D', 'PMinres 2D', 'Chol 3D', 'Minres 3D', 'PMinres 3D'});
	set(gca,'xticklabel', {'jdqr', 'build in', 'PLU', 'PChol', 'PMinres'});
	ylabel('time [s]');
	grid on
	set(gca,'xgrid','off')
end

end

