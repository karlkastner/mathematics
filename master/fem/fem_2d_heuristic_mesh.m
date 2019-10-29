% Thu Jun  7 12:02:45 MSK 2012
% Karl Kästner, Berlin

% Heuristic mesh setup
function [P T Bc Nm G] = fem_2d_heuristic_mesh(L0, x0, tol, order, func)
	% generate an an initial uniform mesh
	n = [3 3];
	[P T Bc] = mesh_2d_uniform(n, L0);
	G = zeros(size(T,1),1);

	% inverse mapping : boundary -> triangle
	[Nm] = fem_2d_element_boundary(P, T, Bc);

	% shift the atom's nucleus
	P = P - ones(size(P,1),1)*x0;

	while (1)
		% estimate the error
		switch (order)
			case {2}
				% get an estimate of the solution
				v_est = feval(func, mesh.P);
				% estimate the error and mark cells for refinement
				[M v_err err_est] = mark_2d_3(P, T, Bc, Nm, v_est);
			case {3}
				mesh = Mesh(P, T, Bc);
				mesh = Promote.promote_2d_3_6(mesh);
				% get an estimate of the solution
				v_est = feval(func,mesh.P);
				% estimate the error and mark cells for refinement
				[M v_err err_est] = mark_2d_6(mesh.P, mesh.T, mesh.Bc, Nm, v_est);
			case {4}
				mesh = Mesh(P, T, Bc);
				mesh = Promote.promote_2d_3_10(mesh);
				% get an estimate of the solution
				v_est = feval(func,mesh.P);
				% estimate the error and mark cells for refinement
				[M v_err err_est] = mark_2d_10(mesh.P, mesh.T, mesh.Bc, Nm, v_est);
			case {5}
				mesh = Mesh(P, T, Bc);
				mesh = Promote.promote_2d_3_15(mesh);
				% get an estimate of the solution
				v_est = feval(func,mesh.P);
				% estimate the error and mark cells for refinement
				[M v_err err_est] = mark_2d_15(mesh.P, mesh.T, mesh.Bc, Nm, v_est);
			case {6}
				mesh = Mesh(P, T, Bc);
				mesh = Promote.promote_2d_3_21(mesh);
				% get an estimate of the solution
				v_est = feval(func,mesh.P);
				% estimate the error and mark cells for refinement
				[area l_boundary h_side s_angle C] = regularity_2d(P,T,Bc);
				[M v_err err_est] = mark_2d_21(mesh.P, mesh.T, mesh.Bc, Nm, v_est, area, h_side, s_angle, C);
			otherwise
				'error'
		end
		if (err_est < tol)
			break;
		end
		% refine triangles
		[P T Bc Nm G] = refine_2d(P, T, Bc, M, Nm, G);
	end % while 1
end % function fem_2d_heuristic_mesh

