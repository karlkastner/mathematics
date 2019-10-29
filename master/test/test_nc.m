% Mon Aug 13 19:43:26 MSK 2012
% Karl Kästner, Berlin

function test_nc()
disturb = 1e-12;
nf = 12;

	e_true = -3*pi^2;
	N_ = [3 4 6 8] + 1; % 6 quad pr.

	bcflag = 1;	
	E = [];
	for idx=1:length(N_)
	n=N_(idx);
	[P T Bc] = mesh_3d_uniform([n n n], [1 1 1]);
	mesh = Mesh_3d(P,T,Bc);
	
	% linear polynomials
	tic();
	f_int = @int_3d_nc_4;
	A = assemble_3d_dphi_dphi_java(mesh, [], f_int);
	B = assemble_3d_phi_phi_java(mesh, [], f_int);
	if (nnz(tril(B,-1)) > 0) error('here'); end
	% apply boundary conditions
	[A B p__] = boundary_3d(A, B, mesh.Bc, bcflag);
	Ta(idx,1) = toc()
	N(idx,1) = size(A,1);
	L = chol(B);
	n_ = min(size(A,1),nf*n)-1;
	LB = -40; UB=0;
	tic()
	E(idx,1) = eigs_lanczos(inv(L')*A*inv(L),n_,n_,1,[],LB,UB);
	%E(idx,1) = eigs_fixed(A,B,1,'SM');
	Tr(idx,1) = toc()

	% quadratic polynomials
	f_int = @int_3d_nc_6;
	tic();
	mesh.promote(2);
	A = assemble_3d_dphi_dphi_java(mesh, [], f_int);
	B = assemble_3d_phi_phi_java(mesh, [], f_int);
	% apply boundary conditions
	[A B p__] = boundary_3d(A, B, mesh.Bc, bcflag);
	Ta(idx,2) = toc()
	N(idx,2) = size(A,1);
	fdx = find(diag(B));
	A=A(fdx,:); A=A(:,fdx);
	B=B(fdx,:); B=B(:,fdx);
	%B = B+disturb*speye(size(B));
	if (nnz(tril(B,-1)) > 0) error('here'); end
	L = chol(B); % + disturb*speye(size(B)));
	LB = -40; UB=0;
	n_ = min(size(A,1),nf*n)-1;
	tic();
	mode = [];
	E(idx,2) = eigs_lanczos(inv(L')*A*inv(L),n_,n_,1,mode,LB,UB);
	Tr(idx,2) = toc()
	%eigs(inv(L')*A*inv(L),[],1,'SM')

	% cubic polynomials
	f_int = @int_3d_nc_8;
	mesh = Mesh_3d(P,T,Bc);
	tic();
	mesh.promote(3);
	Tp(idx,3) = toc()
	tic();
	A = assemble_3d_dphi_dphi_java(mesh, [], f_int);
	B = assemble_3d_phi_phi_java(mesh, [], f_int);
	% apply boundary conditions
	[A B p__] = boundary_3d(A, B, mesh.Bc, bcflag);
	Ta(idx,3) = toc()
	N(idx,3) = size(A,1);
	if (nnz(tril(B,-1)) > 0) error('here'); end
	L = chol(B + disturb*speye(size(B)));
	LB = -40; UB=0;
	n_ = min(size(A,1),nf*n)-1;
	tic();
max(diag(L))
min(diag(L))
	mode = [];
	E(idx,3) = eigs_lanczos(inv(L')*A*inv(L),n_,n_,1,mode,LB,UB);
%	eigs(inv(L')*A*inv(L),[],1,'SM')
	Tr(idx,3) = toc()

	E
	end
	err = E - e_true
	subplot(1,3,1)
	loglog(N,abs(err))
	subplot(1,3,2)
	loglog(N,Tr)
	subplot(1,3,3)
	loglog(N,Ta)
end % function

