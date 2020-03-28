function test_assemble()
%matlabpool('size') > 0
	javaaddpath('.');
	javaaddpath('./fem');
	javaaddpath('/usr/share/java/jama.jar');
%	matlabpool

	T=TestThread(); tic(); T.testThread(2000); toc()    

	n = 300; L0 = [1 1];
	[P T Bc] = mesh_2d_uniform(n, L0);
	P_ = P; T_ = T; Bc_ = Bc;
	int = @int_2d_gauss_3;
	f = Potential_2D_Coulomb; 

	setenv('OMP_NUN_THREADS','1')
	tic()
	A = java_assemble_2d_dphi_dphi(P_, T_, [], int);
	t(1,1) = toc()
	tic();
	V = java_assemble_2d_phi_phi(P_, T_, f, int);
	t(1,2) = toc()

	setenv('OMP_NUN_THREADS','2')
	getenv('OMP_NUN_THREADS')
	tic()
	A_ = java_assemble_2d_dphi_dphi(P_, T_, [], int);
	t(2,1) = toc()
	tic();
	V_ = java_assemble_2d_phi_phi(P_, T_, f, int);
	t(2,2) = toc()
	
	sum(sum((V-V_).^2))
	sum(sum((A-A_).^2))
end % matlabpool	
