
	N_ = [2 4 8 16 32 64 128] % 256 512];

	int = @int_2d_nc_3;
	for idx=1:length(N_)
		n = N_(idx)
		[P T Bc] = mesh_2d_uniform([n n n],[1 1 1]);
		mesh = Mesh_2d(P,T,Bc);
		tic();
		mesh.prefetch();
		Tr(idx,1) = toc();
		tic();
		A = assemble_2d_dphi_dphi_java(mesh, [], int);
		Tr(idx,2) = toc();
		tic();
		B = assemble_2d_phi_phi_java(mesh, [], int);
		Tr(idx,3) = toc();
		N(idx,1:3) = size(A,1);
		Tr
	end

	int = @int_3d_nc_4;
	for idx=1:length(N_)
		n = N_(idx)
		if (n^(3/2) < N_(end))
		[P T Bc] = mesh_3d_uniform([n n n],[1 1 1]);
		mesh = Mesh_3d(P,T,Bc);

		tic();
		mesh.prefetch();
		Tr(idx,4) = toc();
		tic();
		A = assemble_3d_dphi_dphi_java(mesh, [], int);
		Tr(idx,5) = toc();
		tic();
		B = assemble_3d_phi_phi_java(mesh, [], int);
		Tr(idx,6) = toc();
		N(idx,4:6) = size(A,1);
		end
		Tr
	end
	N
	Tr
	loglog(N,Tr,'.-')
	grid on
	set(gca,'minorgrid','none')
