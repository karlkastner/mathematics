% Mon Jun  4 16:38:14 MSK 2012
% Karl Kästner, Berlin

function test_jdqz(mode,n_max)
	path(path,'../jacobi-davidson');
	path(path,'../fem');

	if (nargin<1)
		n_max = 5;
	end
	bcflag = 1;
	%K = [1 10];
	K = [1]; % 10];
	
	% number of grid points
	N = 2.^(2:n_max)'+1;
	
	E_true = -2*pi^2;
	
	% domain size
	L0 = [1 1];

	T_eigs = zeros(size(N,1),1);
	T_jd_man = zeros(size(N,1),1);
	T_jd = zeros(size(N,1),1);

	for kdx=1:length(K)
		k = K(kdx);
	% constant grid loop
	for idx=1:length(N)
		n = N(idx,1)
		if (mode < 2)
			% ordinary eigenvalue problem
			I = speye(n); 
			A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
			A = kron(I,A) + kron(A,I);
			% fake B-matrix for mode 1
			B = kron(I,I);
			k_ = n^2;
		else
			% get the grid
			[P T Bc] = mesh_2d_uniform(n,L0);
		
			% construct 2nd order matrices
			A = assemble_2d_dphi_dphi(P,T,[], 'int_2d_tmp');
			B = assemble_2d_phi_phi(P,T,[], 'int_2d_tmp');
	
			% apply boundary conditions
			[A B P_] = boundary_2d(A, B, Bc, bcflag);
			k_ = size(P_,1);
		end
		if (k_ >= k)
		%k_ = min(size(P_,1), k);
		k_ = k;
		% ARPACK
		tic();
		[v e_eigs] = eigs(A,B,k_,'SM');
		T_eigs(idx,kdx) = toc();

		if (0 == mode)
			% ordinary eigenvalue problem
			% JD - manually implemented
			try
				tic();
%				e_jd  = jacobi_davidson_qr(A, k_, 'SM');
				T_jd_man(idx,kdx) = toc();
			catch cerr
				cerr
				'caught jd_man exception'
				e_jd = zeros(k_,1);
			end
			% JD - Sleijpen
			try
				tic();
				%L = ichol(-A);
				%e_jd  = jdqr(-A, k_, 'SM',L',L);
				%e_jd  = jdqr(A, k_, 'SM');
				[L U] = ilu(A);
				%e_jd  = jdqr(A, k_, 'SM',L,U);
				%T_jd(idx,kdx) = toc();
			catch err	
				'caught jd exception'
			end
		else
			% generalised eigenvalue problem
			% JD
			try
				tic();
				e_jd  = jacobi_davidson_qz(A, B, k_, 'SM');
				T_jd_man(idx,kdx) = toc();
			catch cerr
				
				'caught jdx exception'
				e_jd = zeros(k_,1);
			end
			try
				e_jd  = jdqz(A, B, k_, 'SM');
				T_jd(idx,kdx) = toc();
			catch err	
			end
		end
	end % idx
	end % kdx
	loglog(N,[T_eigs T_jd_man T_jd]);
	%legend('Arpack 1', 'Arpack 10', 'JD 1', 'JD 10');
	legend('Arpack 1', 'JD_man 1', 'JD 1');
	grid on
end
%{
	
	n=40; A=rand(n);
	k=10;
	 I = speye(n);  A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);  A = kron(I,A) + kron(A,I);
	 B = eye(n^2);
	 tic();
	 [v e] = eigs(A,B,k,'SM');% toc(), tic();
	 toc();
	 e = diag(e); 
	%A=A+A'; B=eye(n); tic;
	
	%[e idx] = sort(diag(e));
	%v(:,idx) = v;
	%norm(A*v(:,1))
	%pause
	e'
	%opt.v0=sum(v,2)+0.1*rand(size(v,1),1);
	tic();
	[X_,JORDAN_,Q_,Z_,e_,T_,HISTORY_] = JDQZ(A,B,k,'SM',opt);
	diag(e_)'
	Q_*Z_
	%[X_,JORDAN_,Q_,Z_,e_,T_,HISTORY_] = JDQZ(B,A,k,'LM'); %,opt);
	%1./diag(e_)'
	toc();
	%	[q_ z_ e_(1,idx)] = jacobi_davidson_qz(A,B,v(:,idx)+1e-7*rand(n^2,1), 1);
	pause
	
	tic
	for idx=1:n
		%[q z e] = jacobi_davidson_qz(A,B,sum(v(:,1:10),2), 10); toc
		t=rand(n,1);
		idx
		[q_(:,idx) z_(:,idx) e_(idx)] = jacobi_davidson_qz(A,B,v(:,idx)+0.00*t, 1);
	end
	toc
	q
	q_
	z
	z_
	diag(e)'
	e_
	%, size(q1[q z e] = jacobi_davidson_qz(A,B,sum(v(:,1:10),2), 10); toc, size(q)
%}	
