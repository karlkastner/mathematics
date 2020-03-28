% Tue Jun  5 16:51:06 MSK 2012
% Karl Kästner, Berlin

function test_lanczos(pflag)
	'1D'
	k  = 20;
	N_ = 40*[1 1 1];
	n=N_(1);
	A=(n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	tic(); E = eigs(A,[],n-1); t(1,1)=toc();
	%E_true = zeros(k,1); E=unique(single(E)); E_true(1:length(E)) = E;
	tic(); E_ = eig_symmetric(A,1); t(1,2)=toc();
	%E=sort(E); E_=sort(E_); E_=E_(1:end-1); norm(E-E_)
	E =[E;0];
	[E E_ E_-E] 
	t

	'2D'
	%n=round(n_^(3/2));
	n = N_(2);
	m = 2*2*n;
	A=(n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	I = speye(n);
	A = kron(A,I) + kron(I,A);
	tic(); E=eigs(A,[],k,'SM'); t(2,1)=toc();
	E_true = zeros(k,1); E__=unique(single(E)); E__=sort(E__,1,'descend'); E_true(1:length(E__)) = E__;
	tic(); E_ = eigs_lanczos(A, k, m, 'SM'); t(2,2)=toc();
	[E E_ E_-E_true] 
	t
	'3D'
	n = N_(3);
	m = 2*3*n;
	A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	I = speye(n);
	A = kron(kron(A,I),I) + kron(kron(I,A),I) + kron(kron(I,I),A);
	tic(); E  = eigs(A,[],k,'SM'); t(3,1)=toc();
	E_true = zeros(k,1); E__=unique(single(E)); E__=sort(E__,1,'descend'); E_true(1:length(E__)) = E__;
	tic(); E_ = eigs_lanczos(A, k, m, 'SM'); t(3,2)=toc();
	[E E_ E_-E_true] 
	t

	figure(1); clf; subplot(1.5,1.5,1);
	bar(t);	
	title(['Run Times For Computing ' num2str(k) ' Eigenvalues of a ' num2str(n) '^d sized Poisson System']);
	set(gca,'xticklabel', {'1D', '2D', '3D'});
	legend('location','northwest','Arpack','Lanczos')
	ylabel('time [s]');
	grid on
	set(gca,'xgrid','off')
	if (nargin() > 0 && pflag)
		print -depsc eigensolver-poisson-run-times-ii.eps
	end
end % test_lanczos

% duplicate eigenvalus are not found
% no eigenvectors (can be found by shift and invert afterwards)
% can be sped up by factor of 10 (needs porting to java)
% does it apply for matrices generated from meshes?

