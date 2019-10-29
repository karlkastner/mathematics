% Jul  2 20:17 
% Karl Kästner, Berlin


function test_compare_solvers()
	addpath('../jacobi-davidson')
	addpath('../cullum')

        k = 16;
	
	n = 100*[1 2];
%	[Err t] = compare(1,n)
	[Err t] = compare(4,n)
%	[Err t] = compare(16,n)
%	[Err t] = compare(64,n)

	n = 15*[1 2 4];
%	[Err t] = compare(1,n)
	[Err t] = compare(4,n)
%	[Err t] = compare(16,n)
%	[Err t] = compare(64,n)

end

function [Err t] = compare(k,n)
	A = poisson(n); %/(10*max(n)^2);

	tic();
	E1 = eigs(A,[], k,'SM'); %,opts);
	t(1) = toc();
	Err(1) = 0;


	tic();
	jdopts.maxit = size(A,1);
%	opt.jmin = 1;
%	opt.LS_Tol = 1e-2;
%	opt.jmax = ;
	E2 = jdqr(A, k, 'SM', jdopts);
	t(2) = toc();
%	Err(2) = norm(E1-E2);
	Err(2) = 0;

	tic();
	UB = 0;
	%LB = -20*k^2;
	LB = 2*min(E1);
%	KMAX = min(k*max(n),prod(n));
%	KMAX = min(round(k*max(n)*length(n)/2),prod(n));
	KMAX = min(round(k*sum(n)),prod(n));
	NMEV = KMAX;
	E3 = eigs_lanczos_wrapper(A, KMAX, NMEV, LB, UB);
	t(3) = toc();
	E3 = sort(E3,'descend');
	l = min(length(E1),length(E3))
	Err(3) = norm((E1(1:l) - E3(1:l)));
%	[E1(1:l)  E3(1:l)]
end
