% Fri Jul 27 15:42:39 MSK 2012
% Karl Kästner, Berlin

	
function test_Mesh_3d()
	format compact;
	 %opengl neverselect;

	plotflag = 1;
	i_max = 4;
	p = 0;
	m = 1;

	n =  [1 1 1]*2;
	L0 = [1 1 1];
	x0 = [0 0 0];

	% initial mesh
	[P T Bc] = mesh_3d_uniform(n,L0,x0);
	[v_sum0 a_sum0 h_eff_max volume area h_eff] = regularity_3d(P,T,Bc);
	clf();
	% loading mesh into the tree
	tree = Tree_3d(P, T, Bc);

	% generate the mesh
	mesh = tree.generate_mesh();

	% check, wether the mesh was successfully loaded
	[P_ T_ Bc_] = get_mesh_arrays(mesh);
	[v_sum a_sum h_eff_max volume area h_eff] = regularity_3d(P_, T_, Bc_);
	T = sortrows(sort(T,2));
	T_ = sortrows(sort(T_,2));
	Bc = sortrows(sort(Bc,2));
	Bc_ = sortrows(sort(Bc_,2));
	tree.check_neighbour();
	[norm(P-P_) norm(T-T_) norm(Bc-Bc_) v_sum-v_sum0 a_sum-a_sum0]
%	clf();
%	subplot(2,2,1)
%	tetramesh(T_,P_);

	% refining
	M_ = [2 16];
	for wdx=1:100
		wdx
		clf();
		[P T Bc] = mesh_3d_uniform(n,L0,x0);
		tree = Tree_3d(P,T,Bc);
		for idx=1:i_max
			N(idx,1) = size(P,1);
			%n_ = sum(abs([P(T(:,1),:) P(T(:,2),:) P(T(:,3),:) P(T(:,4),:)]),2); [min_ imin] = min(n_); M = imin;
			%M = randi(size(T,1),1,randi(round(size(T,1)/4)))
			%M=randi(size(T,1)) % randi(size(T,1))]
			%if (idx >= i_max-1) M=randi(size(T,1),1,randi(size(T,1))); end
			if (idx > 1) M =M_(idx-1); tree.refine(M); end
			mesh = tree.generate_mesh();
			mesh.element_neighbours();
			[P T Bc N] = get_mesh_arrays(mesh);
			[v_sum a_sum h_eff_max volume area h_eff] = regularity_3d(P_, T_, Bc_);
			[v_sum-v_sum0 a_sum-a_sum0]
			[T sdx] = sort(T,2);
%			[ (1:length(P))' P]
			l = size(unique(T,'rows'))
			size(T,1)
			l = size(unique(P,'rows'))
			size(P,1)
			%N(sdx,:) = N;
			%[T sdx] = sortrows(T); %sort(T,2));
			%N = sort(N,2); 
			%N=N(sdx,:);
%			[T N]

			if (idx > 2)
			if (plotflag)
				figure(1);
				subplot(sqrt(i_max),sqrt(i_max),idx);
				cla();
				display_3d(P,T,Bc,5+8+p+16);
				axis([0 1 0 1 0 1]-x0(1))
			
				figure(2);
				subplot(sqrt(i_max),sqrt(i_max),idx);
				cla;
				size(Bc,1)
				display_3d(P,T,Bc,6+8);
				axis([0 1 0 1 0 1]-x0(1))
				pause
			end
			end
		end % for idx
	end % for wdx
%	figure(3)
%	loglog(N,[Tr H_min])
	% TODO, 
%	legend(L);	
	% refine radomly (TODO: no closure by now)
%	M = randi(size(T,1),round(size(T,1)/8));
%	M = fliplr((1:size(T,1)));
%	tree.refine(M)
%	mesh = tree.generate_mesh();
%	P_ = double(mesh.P); P_ = P_(1:mesh.np,:);
%	T_ = double(mesh.T); T_ = T_(1:mesh.nt,:);
%	Bc_ = double(mesh.Bc); Bc_ = Bc_(1:mesh.nb,:);
%	[v_sum a_sum h_eff_max volume area h_eff] = regularity_3d(P_,T_,Bc_);
%	v_sum
%	a_sum

end % test_Mesh_3d

