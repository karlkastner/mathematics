% Mon Aug 13 03:17:13 MSK 2012
% Karl Kästner, Berlin

	n = 20;
	N = [n 2*n 4*n];
	A = poisson(N);

	size(A,1)
	tic()
	[v e] = eigs(A,[],1,'SM');
	e
	toc()
	tic()
	%eigs(A,[],1,'LA')
	eigs(-A,[],1,'SA')
	toc()
	KMAX=2*n*4;
	NMEV=KMAX;
	LB=-30; UB=0;
	tic();
	E = eigs_lanczos(A, KMAX, NMEV, 1, 'SA', LB, UB)
	toc()
	tic()
	E = eigs_lanczos(A, KMAX, NMEV, 1, 'LA', LB, UB)
	toc()
	tic()
	I = speye(size(A));
	b = rand(size(A,1),1); b=b/norm(b);
	tol = sqrt(eps);
	MAXIT = 8*round(size(A,1).^(1/3))
	x=minres(A-E(1)*I,b,tol,MAXIT);
	x = x/norm(x);
	[norm(A*x-E(1)*I*x) norm(x-v)]
	toc()
%{
	IDO = 0;
	BMAT = 'I';
	N = size(A,1);
	WHICH = 'LA';
	NEV = 1;
	TOL = sqrt(eps);
	RESID = rand(N,1);
	NCV = 5;
	V = zeros(N,NCV);
	LDV = N;
	IPARAM = zeros(11,1);
	IPTR = zeros(15,1);
	WORKD = zeros(N,3);

	INFO = 0;

	arpackc('dsaupd',IDO,BMAT,N,WHICH,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,WORKL,LWORKL,INFO);
	% IDO: reverse communication parameter, initialized to 0. {int}
	%    BMAT: 'I' for standard problem, 'G' for generalized. {char}
	%    N: size of problem. {int}
   	% WHICH: 'LM','SM','LA','SA','BE'. {length 2 char}
	%    NEV: number of eigenvalues requested. {int}
	%    TOL: convergence tolerance. Default is eps/2. {double}
	%    RESID: may be initialized to start vector. {length N double}
	%    NCV: number of Lanczos vectors. {int}
	%    V: Lanczos basis vectors. {length N*NCV double}
	%    LDV: leading dimension of V. {int}
	%    IPARAM: {length 11 int}
	%    IPNTR: {length 15 int}
	%    WORKD: {length 3*N double}
    WORKL: {length LWORKL double}
    LWORKL: length of WORKL >= NCV^2+8*NCV {int}
    INFO {int}
%}

