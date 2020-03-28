% Wed Jun 13 11:57:07 MSK 2012
% Karl Kästner, Berlin

function test_2d_refine_structural()
	k = 12;
	order = 2;
	% initial grid
	[P T B X] = mesh_2d_uniform([3 5], [160 160]);
	X1 = X{1} - 80;
	X2 = X{2} - 80;
	X = {X1, X2};
	figure(1); clf
	figure(2); clf
	for idx=1:k
		% set up structural grid
		[P T Bc X] = mesh_2d_uniform([], X);
		% inverse mapping : boundary -> triangle
		Nm = fem_2d_element_boundary(P, T, Bc);
		% get estimated (anlytic) solution
		v = hydrogen_2d_analytic(P);
		% estimate the error and mark cells for refinement
		[M v_err err_est] = fem_2d_mark_3(P, T, Bc, Nm, v);
		% structural refinement
		X = fem_2d_refine_structural(X,v_err,order);
		% display solution and grid
		figure(1);
		subplot(ceil(sqrt(k)), ceil(sqrt(k)), idx);
		%display_2d(P, T, Bc, 0, [v(T(:,1),1).^2 v(T(:,2),1).^2 v(T(:,3),1).^2], [], 'EdgeColor', 'none');
		display_2d(P,T,Bc,0,log10(v_err),[],'EdgeColor','none');
		axis equal, axis tight
		figure(2);
		subplot(ceil(sqrt(k)), ceil(sqrt(k)), idx);
		display_2d(P,T,Bc,8);
		axis equal, axis tight
		drawnow
	end % for idx
end % function

