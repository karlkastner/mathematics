% Wed Aug 15 22:45:15 MSK 2012
% Karl Kästner, Berlin


function test_2d_eigensolver()

	nf = 1;
	n_jd=2;
	folder = '../dat/2d-mesh-series/';
	t_max = 30;

	bcflag = 1;
	name_c = ls(folder,'-1');
	name_c = regexp(name_c,'\n','split')
	name_c = name_c(2:end);

	lw = 2;
	ms = 12;
	shift = -3;
	eig_mode = 'SM';

	k = 1;
	int = @int_2d_nc_3;

	Tr = zeros(0,0);
	N  = zeros(0,1);
	E  = zeros(0,0);

%	arpack_opt.p = 5;
	T_max = zeros(1,3);

	for idx=1:length(name_c)
		if (11==idx) continue; end
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
			potential = Potential_2D_Coulomb;
			% load mesh
			mesh = Mesh_2d(s.P,s.T,s.Bc);
			% assemble the discretisation matrices
			A = assemble_2d_dphi_dphi_java(mesh,[],int);
			B = assemble_2d_phi_phi_java(mesh,[],int);
			V = assemble_2d_phi_phi_java(mesh,potential,int);
			A = -0.5*A + V;
			% apply the boundary conditions
			[A B] = boundary_2d(A,B,mesh.Bc,bcflag);

			% resave
			s.A = A;
			s.B = B;
		
			% resave
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
		
		l = 1;
		% 1 solve with eigs (manually shifted 'SM')
		if (T_max(l) < t_max)
			'eigs'
			tic();
			% todo, try also SM
			e = eigs(A-shift*B,B,k,eig_mode, arpack_opt);
			E(idx,l) = e + shift;
			Tr(idx,l) = toc();
			T_max(l) = max(T_max(l),Tr(idx,l));
		end
		[log10(N) Tr E]
%{
		% 2 solve with eigs (internal shift)
		l=l+1;
		if (1 == idx || Tr(idx-1,l) < t_max)
			tic();
			% todo, try also SM
			e = eigs(A,B,k,shift,arpack_opt);
			E(idx,l) = e;
			Tr(idx,l) = toc();
		end
		[log10(N) Tr E]

		% 4 solve with eigs (transformed to standard problem)
		l=l+1;
		if (1 == idx || Tr(idx-1,l) < t_max)
			tic();
			[v e flag] = eigs(A_,[],k,eig_mode,arpack_opt);
			E(idx,l) = e+shift;
			Tr(idx,l) = toc();
		end
		[log10(N) Tr E]

% {
		% 3 solve with eigs (smallest algebraic, manual shift)
		l=l+1;
		if (1 == idx || Tr(idx-1,l) < t_max)
			tic();
			% todo, try also SM
			e = eigs(A-shift*B,B,k,'SA',arpack_opt);
			E(idx,l) = e+shift;
			Tr(idx,l) = toc();
		end
		[log10(N) Tr E]

		% 5 solve with eigs (transformed to standard problem, SA)
		l=l+1;
		if (1 == idx || Tr(idx-1,l) < t_max)
			tic();
			% todo, try also SM
			e = eigs(A_,[],k,'SA',arpack_opt);
			E(idx,l) = e+shift;
			Tr(idx,l) = toc();
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
%				jdots.jmin = 40;
%				jdots.jmax = 20;
				tic();
%			       	[v e] = jdqz(A-shift*B, B, k, eig_mode, jdopts);
				e = NaN;
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
%				e = eigs_lanczos(A_, n, n, k, 'SA', -10, 10);
%				e = eigs_lanczos(A_, n, n, k, 'SA', [], []);
%				e = eigs_lanczos(A_, n, n, k, eig_mode, [], []);
				e = e+shift;
				Tr(idx,l) = toc();
				T_max(l) = max(T_max(l),Tr(idx,l));
				if (isempty(e) || isnan(e) || abs(e - E(idx,1)) > 1e-3)
					nf = 2*nf;
					if (n >= 10*size(A,1) || Tr(idx,l) > t_max)
						E(idx,l) = NaN;
						break;
					end
				else
					E(idx,l) = e;
					break;
				end
			end
			% todo, time minres to find the corresponding eigenvector (with ilu)
		end
		[log10(N) Tr E]
		Ap = poisson(round(size(A,1)^(1/3))*[1 1 1]);
		log10([condest(A) condest(A_) condest(Ap)])
%		[v e] = eig(full(A),full(B)); cond(v)
%		[v e] = eig(full(A_)); cond(v)
	end	
	
	[log10(N) Tr]
	[log10(N) E]
	loglog(N,Tr,'linewidth',lw,'markersize',ms)
	preparePrint();
	print -depsc ../img/2d-benchmark-eigensolver-hydrogen.eps
	legend('Arpack','Jacobi-Davidson','Lanczos')
	xlabel('n : number of grid points')
	ylabel('run time [s]')

end

