% Wed Jun 13 21:45:12 MSK 2012
% Karl Kästner, Berlin

function test_refine_2d()
	% initial mesh
	n = [2 2]; L = [1 1];
	[P T Bc] = mesh_2d_uniform(n,L);
	for idx=1:5
		[P T Bc] = mesh_2d_uniform(n, L);
		% check domain volume and surface area and degeneracy
		[a_sum l_sum h_eff_max] = regularity_2d(P,T,Bc);
		[a_sum l_sum h_eff_max]
		length(P)
		% uniform refinement
		% TODO random refinement
		n = n*2;
	end
end

