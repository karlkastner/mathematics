% Sun Aug 12 14:32:18 MSK 2012
% Jul 11 18:54
% Karl Kästner, Berlin

function test_mg()
	tol = 1e-10;
	n_max = 5;
	int = @int_3d_gauss_4;

	% initial mesh
	L = 1;
	n = 4;
	[P T Bc] = mesh_3d_uniform([n n n], [L L L]);
	tree = Tree_3d(P,T,Bc);

	for idx=1:n_max
		mesh = tree.generate_mesh();
		% assemble the B-matrix
		B = assemble_3d_phi_phi_java(mesh,[],int);
		% apply boundary conditions
		[void B] = boundary_3d(sparse(size(B)), B, Bc);

%{
		% assemble the projection matrices
		if (idx > 1)
			M = mg_mat(P,Q,poly);
			M_cell{idx} = M; % what about back projection ?
		else
			M_cell{idx} = [];
		end
		B_cell{idx} = B;
%}

		% solve the system
		n = size(B,1);
		N(idx,1) = n;
		b = randn(n,1);
		%b = zeros(n,1); b(round(n/2)) = 1;
		tic();
		x_mg = 0; %mg_solve(B_cell, M_cell, b, b, tol);
		Tr(idx,1) = toc();
		tic();
		x_jacobi  = jacobi(B, b, b, tol);
		Tr(idx,2) = toc();
		tic();
		x_gs  = gauss_seidel(B, b, b, tol);
		Tr(idx,3) = toc();
		tic();
		tic();
		x_mr  = minres(B, b, tol, size(B,1), [], [], b);
		Tr(idx,4) = toc();
		tic();
		p = symamd(B);
		x(p,1) = B(p,p) \ b(p);
		Tr(idx,5) = toc();
		norm(x-x_mg)
		norm(x-x_jacobi)
		norm(x-x_gs)
		norm(x-x_mr)
		Tr

		% uniformly refine the mesth
		if (idx < n_max)
			[P T Bc] = get_mesh_arrays(mesh);
			R = (1:size(T,1));
			tree.refine(R);
		end
	end
end % function test_mg

function [x flag nrvec r] = gauss_seidel(A, b, x, tol)
	p = 0.0;
	U = triu(A,+1);
	L = tril(A);
	maxiter = 10*size(A,1);
	nrvec = [];
	idx=1;
	while (1)
		r = A*x - b;
		nrs = r'*r;
		nrvec(idx) = nrs; 
		if (nrs < tol^2)
			flag = 0;
			break;
		end
		x = p*x + (1-p)* L \ (b - U*x);
		idx = idx+1;
		if (idx > maxiter)
			flag = 1;
			break;
		end
	end
end

function [x flag nrvec r] = jacobi(A, b, x, tol)
	p = 0.2;
	D = diag(A);
	iD = 1./D;
	R = A - diag(sparse(D));
	maxiter = 10*size(A,1);
	nrvec = [];
	idx=1;
	while (1)
		r = A*x - b;
		nrs = r'*r;
		nrvec(idx) = nrs; 
		if (nrs < tol^2)
			flag = 0;
			break;
		end
		x = p*x + (1-p)*iD.*(b - R*x);
		idx = idx+1;
		if (idx > maxiter)
			flag = 1;
			break;
		end
	end
end

function x = v_step(Ac, Mc, b, x, idx)
	% one jacobi iteration smoothing step
	A = Ac{idx}
	D = diag(A);
	R = A - D;
	iD = 1./D;
	r = A*x-b;
	x = iD.*(b - R*x);
	idx=idx+1;
	% check for terminal condition
	if (idx > 1)
		% contract
		M = Mc(idx);
		x_ = M*x;
		b_ = M*b;
		[x_] = v_step(Ac, Mc, b_, x_, idx-1);
		% expand
		x = M'*x_;
		% apply changes
		x = iD.*(b - R*x);
	end
end % function

