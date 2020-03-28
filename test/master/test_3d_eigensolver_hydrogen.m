% Wed Aug 15 23:14:31 MSK 2012
% Karl Kästner, Berlin


function test_3d_eigensolver_hydrogen()

	nf = 1;
	n_jd=2;
	folder = '../dat/3d-mesh-series/';
	t_max = 30;

	bcflag = 1;
	name_c = ls(folder,'-1');
	name_c = regexp(name_c,'\n','split')
	name_c = name_c(2:end);

	lw = 2;
	ms = 12;
	shift = -0.55;
	eig_mode = 'SM';

	k = 1;
	int = @int_3d_nc_4;

	Tr = zeros(0,0);
	N  = zeros(0,1);
	E  = zeros(0,0);

	%arpack_opt.p = 5;
	T_max = zeros(1,3);

	for idx=1:length(name_c)
		% load data
		try
			filename = [folder '/' name_c{idx}]
			s = load(filename);
		catch err
			continue
		end
		
		% check if the discretisation matrices where saved
		if (~isfield(s,'A'))
			'setting up the discretisation matrices'
			potential = Potential_3D_Coulomb;
			% load mesh
			mesh = Mesh_3d(s.P,s.T,s.Bc);
			% assemble the discretisation matrices
			A = assemble_3d_dphi_dphi_java(mesh,[],int);
			B = assemble_3d_phi_phi_java(mesh,[],int);
			V = assemble_3d_phi_phi_java(mesh,potential,int);
			A = -0.5*A + V;
			% apply the boundary conditions
			[A B] = boundary_3d(A,B,s.Bc,bcflag);
		
			% resave
			s.A = A;
			s.B = B;
			save(filename, '-struct', 's');
		else
			'loading the discretisation matrices'
			% load variables
			A = s.A;
			B = s.B;
		end

		L = chol(B); iL = inv(L);
		A_ = iL'*(A-shift*B)*iL;
		A_ = 0.5*(A_+A_');

		% numebr of grid points
		n = size(A,1);
		N(idx,1) = n;
		arpack_opt.maxit = n;
		
		% solve with eigs
		l = 1;
		if (T_max(l) < t_max)
			'eigs'
			tic();
			% todo, try also SM
			e = eigs(A-shift*B,B,k,eig_mode, arpack_opt);
			E(idx,l) = e+shift;
			Tr(idx,l) = toc();
			T_max(l) = max(T_max(l),Tr(idx,l));
		end

%{
		% solve with eigs
		if (1 == idx || Tr(idx-1,2) < t_max)
			tic();
			% todo, try also SM
			E(idx,2) = eigs(A_,[],k,'SA');
			Tr(idx,2) = toc();
		end
		[log10(N) Tr E]
%}

		% solve with jacobi davidson
		l = l+1;
		if (T_max(l) < t_max)
		'jd'
			while (1)
				jdopts.maxiter  = n_jd;
				jdopts.LS_MaxIt = n_jd;
				jdopts.tol = sqrt(eps);
%				jdots.jmax = 10;
%				jdopts.LS_MaxIt = max(10,size(A,1));
				tic();
				[v e] = jdqz(A-shift*B, B, k, eig_mode, jdopts);
				Tr(idx,l) = toc();
				T_max(l) = max(T_max(l),Tr(idx,l));
				e =  real(e)+shift;
				if (isempty(e) || isnan(e) || abs(e - E(idx,1)) > 1e-3)
					n_jd = 2*n_jd;
					if (n_jd >= size(A,1) || T_max(l) > t_max)
						E(idx,l) = NaN;
						break;
					end
				else
					E(idx,l) = e;
					break;
				end
%			       	[v e] = jdqz(A, B, k, 'SR', jdopts);
%				if (isempty(e)) e = NaN; end
%				E(idx,l) = real(e);
			end
		end
		[log10(N) Tr E]

		% solve with lanczos
		l=l+1;
		if (T_max(l) < t_max)
		'lanczos'
			while (1)
				n = size(A,1);
				n = min(n, round(nf*n^(1/2)));
				tic();
				e = eigs_lanczos(A_, n, n, k, 'SA', -1-shift, 0-shift);
%				e = eigs_lanczos(A_, n, n, k, 'SA', [], []);
%				e = eigs_lanczos(A_, n, n, k, eig_mode, -2, 0);
				e = e+shift;
				Tr(idx,l) = toc();
				T_max(l) = max(T_max(l),Tr(idx,l));
				if (isempty(e) || isnan(e) || abs(e - E(idx,1)) > 1e-3)
					nf = 2*nf;
					if (n >= 2*size(A,1) || T_max(l) > t_max)
						E(idx,l) = NaN;
						break;
					end
				else
					E(idx,l) = e;
					break;
				end
			end % while 1
			% todo, time minres to find the corresponding eigenvector (with ilu)
		end
		[log10(N) Tr E]
%		Ap = poisson(round(size(A,1)^(1/3))*[1 1 1]);
%		log10([condest(A) condest(A_) condest(Ap)])
%		[v e] = eig(full(A),full(B)); cond(v)
%		[v e] = eig(full(A_)); cond(v)
	end	% for idx
	
	[log10(N) Tr]
	[log10(N) E]
	loglog(N,Tr,'linewidth',lw,'markersize',ms)
	preparePrint();
	print -depsc ../img/3d-benchmark-eigensolver-hydrogen.eps
	legend('Arpack','Jacobi-Davidson','Lanczos')
	xlabel('n : number of grid points')
	ylabel('run time [s]')

end

