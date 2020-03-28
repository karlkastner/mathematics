% Tue May 15 14:23:51 MSK 2012
% Karl Kästner, Berlin

function test_fem_2d_adaptive()
	path(path,'fem')
	n = 3;
	L0 = 1;
	E_true = -2*pi^2;

	% get the initial grid
	[P T Bc] = mesh_2d_uniform(n,L0);
	% get initial neighbourhood relations
	[Nm] = fem_2d_element_boundary(P, T, Bc);
	G=zeros(size(T,1),1);

	P3 = P; T3 = T; Bc3=Bc; G3 = G; Nm3 = Nm;

	for idx=1:7
		N(idx,1) = size(P3,1)
		N(idx,2) = size(P,1)

		% second order

		% construct discritisation matrices
		A = assemble_2d_3_dphi_dphi(P3, T3);
		B = assemble_2d_3_phi_phi(P3, T3, 'int_2d_gauss');
		% apply boundary conditions
		[A B] = boundary_2d(A, B, Bc3);

		[v e] = eigs(A,B,1,'SM');
		E(idx,1) = e;
		% mark cells for refinement
		[M err err_max dV C] = fem_2d_mark_3(P3, T3, Bc3, Nm3, v);
		Err_est(idx,1) = err_max;
		% refine the mesh
		[P3 T3 Bc3 Nm3 G3] = fem_2d_refine(P3, T3, Bc3, M, Nm3, G3);

		% third order

		% promote the triangles to 6-point
		[P6 T6 Bc6] = promote_2d_3_6(P, T, Bc);
	
		% construct 3rd order matrices
		A = assemble_2d_6_dphi_dphi_gauss(P6,T6);
		B = assemble_2d_6_phi_phi_gauss(P6,T6);
	
		% apply boundary conditions
		[A B] = boundary_2d(A, B, Bc6);
	
		% solve the 3rd order system
		[v e] = eigs(A,B,1,'SM');
		E(idx,2) = e
	
		% mark cells for refinement
		[M err err_max dV C] = fem_2d_mark_6(P6, T6, Bc6, Nm, v);
		Err_est(idx,2) = err_max;
	
		% refine the mesh
		[P T Bc Nm G] = fem_2d_refine(P, T, Bc, M, Nm, G);

		(E - E_true)/E_true
		Err_est
	end % for idx

	Err = abs((E - E_true)/E_true);
	disp([log(1./N) Err Err_est])
	% plot
	loglog(N(:,1),Err(:,1),'.-b'); hold on
	loglog(N(:,1),Err_est(:,1),'.-g');
	loglog(N(:,2),Err(:,2),'.-r');
	loglog(N(:,2),Err_est(:,2),'.-m'); hold off
	xlabel('N : total number of mesh points');
	ylabel('absolute error of first eigenvalue');
end % test_fem_2d_adaptive

