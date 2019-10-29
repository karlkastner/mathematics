% Wed Jun 27 16:13:30 MSK 2012
% Karl Kästner, Berlin

% compares the Lanczos eigensolver with the Arpack eigensolver
% at the discrete Laplacian operator
function test_eigs_lanczos()

%	n = 1000;
%       n = 400*[1 2];
	n = 12*[1 2 4];

	% generate a poisson matrix
	A = poisson(n);
	% number of eigenvalues to compute
        k = 10;

        tic()
	% call ARPACK
        E1 = eigs(A,[],k,'SM');
        t(1) = toc();

        tic()
	% bounds of the bisection routine
	LB = min(E1)-1;
	UB = max(E1)+1;
	% matrix rank
        N = size(A,1);
	% size of the T-matrix
	KMAX = min(N,k*max(n));
	% eigenvalues of the T-matrix to be computed
	NMEV = KMAX-1;
	% call FORTRAN eigensolver
	E2 = eigs_lanczos(A, KMAX, NMEV, LB, UB);
        t(2) = toc();

	% display run times
	fprintf('run time eigs   %f\n',t(1));
	fprintf('runtime Lanczos %f\n',t(2));
	fprintf('\n');
	
	% display eigenvalues
	E1 = sort(E1,'descend');
	E2 = sort(E2,'descend');
	k2 = min(length(E1), length(E2));
	fprintf(1,'number of desired eigenvalues:   %d\n',k);
	fprintf(1,'number of converged eigenvalues: %d\n',k2);
	fprintf('\n');
	fprintf('E_eigs E_lanczos difference\n');
	fprintf(1,'%+1.4e %+1.4e %+1.4e\n',[E1(1:k2) E2(1:k2) E1(1:k2)-E2(1:k2)]');
end

