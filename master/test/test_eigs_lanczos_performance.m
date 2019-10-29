% Mon Aug 13 19:47:03 MSK 2012
% Karl Kästner, Berlin

function test_eigs_lanczos_performance()
	nf = 12;
	LB = -40; UB=0; mode = [];
	N_ = [2 4 8 16]; % 16 32];
	bcflag = 1;
	e_true = -3*pi^2;
	
	for idx=1:length(N_)
		n = N_(idx)
		A = poisson([n 2*n 4*n]);
		N(idx,1) = size(A,1)
		tic()
		E(idx,1) = eigs_lanczos(A, nf*n, nf*n, 1, 'LA', [], []);
		Tr(idx,1) = toc()
		
		tic();
		[P T Bc] = mesh_3d_uniform([2*n 2*n 2*n], [1 1 1]);
		mesh = Mesh_3d(P,T,Bc);
		f_int = @int_3d_nc_4;
		A = assemble_3d_dphi_dphi_java(mesh, [], f_int);
		B = assemble_3d_phi_phi_java(mesh, [], f_int);
		% apply boundary conditions
		[A B p__] = boundary_3d(A, B, mesh.Bc, bcflag);
		Ta(idx,2) = toc();
		N(idx,2) = size(A,1);
		L = chol(B);
		n_ = min(size(A,1),nf*n)-1;
		tic();
		E(idx,1) = eigs_lanczos(inv(L')*A*inv(L),n_,n_,1,mode,LB,UB);
		Tr(idx,2) = toc();
	end
	subplot(1,2,1)
	loglog([N N(:,2)],[Tr Ta(:,2)],'.-')
	subplot(1,2,2)
	Err = E - e_true;
	loglog(N,abs(Err),'.-')

end % test_eigs_lanczos_performance()
	
