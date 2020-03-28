% Wed Jun 13 21:14:27 MSK 2012
% Karl Kästner, Berlin

% heaxhedral

function test_refine_3d()
	% initial mesh
	n = 3; L = 1;
	e_true = -3*pi^2;
	[P T Bc] = mesh_3d_uniform(n,L);
	tree = Tree_3d(P,T,Bc);
	for idx=1:4
		%[P T Bc] = mesh_3d_uniform(n,L);
		% check domain volume and surface area and degeneracy
		%[v_sum a_sum h_eff_max] = regularity_3d(P,T,Bc);
		%[v_sum a_sum h_eff_max]
		%length(P)
		% uniform refinement
		%[P T Bc] = refine_3d(P, T, Bc, 1:size(T,1));
		% TODO random refinement
		%n = n*2;
		mesh = tree.generate_mesh();
%		mesh.promote_4_10();
		int = @int_3d_gauss_11;
		A = assemble_3d_dphi_dphi_java(mesh,[],int);
		B = assemble_3d_phi_phi_java(mesh,[],int);
		[P T Bc] = get_mesh_arrays(mesh);
		[A B p__] = boundary_3d(A,B,Bc,1);
		if (1 == size(A,1)) A=full(A); B=full(B); end
		e = eigs(A,B,1,'SM');
		err = e - e_true;
		mesh.prefetch();
		[idx mesh.np size(A,1)]
		[ min(mesh.h_max) max(mesh.h_max) e err]

		M = (1:mesh.nt)';
		tree.refine(M);
		[v_sum a_sum h_eff_max volume area h_eff h_max] = regularity_3d(P,T,Bc);
		v_sum
		a_sum
		
%		subplot(2,2,idx)
		[P T Bc] = get_mesh_arrays(mesh);
%		display_3d(P,T,Bc,1);
	end
end


