function test_inverse_problem()
		bcflag = 1;

		L0 = 2*[ 5  5]
		n  = 1000*[1 1]

		for jdx=2:10;
		n = (2^jdx + 1)*[1 1];

		f_int = @int_2d_gauss_1;
%		f_pot = @Potential_2D_Coulomb
		f_pot = javaObject('Potential_2D_Coulomb');
		% get the grid
		[P T Bc] = mesh_2d_uniform(n,L0);
		mesh = javaObject('Mesh_2d', P, T, Bc);
		mesh.prefetch();
		mesh.element_neighbours();

		A = assemble_2d_dphi_dphi_java(mesh, [], f_int);
		B = assemble_2d_phi_phi_java(mesh, [], f_int);
		V = assemble_2d_phi_phi_java(mesh, f_pot, f_int);
		A = -0.5*A+V;

		Bc = mesh.Bc;
		[A B p] = boundary_2d(A, B, Bc, bcflag);

		mu = -2

%		tic()
%		e = eigs(A-mu*B,B,1,'LM');
%		e + mu
%		toc()

		% test definiteness of A-muB
%		eigs(A-mu*B,[],1,'SA')
%		eigs(A-mu*B,[],1,'LA')
%		eigs(B,[],1,'SA')
%		eigs(B,[],1,'LA')
%		eps*condest(A-mu*B)
%		eps*condest(B)

		tic()
		e = eigs(A-mu*B,B,1,'SM');
		e + mu
		T_(jdx,1) = toc()

%		tic()
%		e = eigs(B,A-mu*B,1,'LM');
%		1/e + mu
%		toc()


%		opts.p = 20;
%		opts.v0 = ones(size(A,2),1);
%		opts.tol = 1e-7;
		opts.isreal = 1;
		opts.issym = 1;
		opts.tol = 1e-4;
%		opts

		tic()
		nc=0;
		mc=0;
		AmuB = A-mu*B;
%		L = ichol(AmuB);
		popts.michol='on';
%		popts.type = 'ict';
%		popts.droptol = 0.125;
		L = ichol(AmuB,popts);
%		L = speye(size(B));
		[v e flag] = eigs(@(x) Afunc(AmuB,B,x,L),size(A,2),1,'SM',opts);
		e + mu
		nc
		mc
		flag
		T_(jdx,2:4) = [toc() mc nc]

% multigrid, AMG, sparse approximate inverse

% ichol precontioner (ict droptols do not improve)
%   1.0e+03 *
%
%         0         0         0         0
   % 0.0000    0.0000    0.0270    0.0100
  %  0.0000    0.0000    0.0420    0.0210
 %   0.0000    0.0001    0.0840    0.0210
%    0.0001    0.0001    0.1460    0.0210
   % 0.0002    0.0003    0.2800    0.0210
  %  0.0008    0.0016    0.5550    0.0210
 %   0.0038    0.0118    1.1130    0.0210
%    0.0154    0.1062    2.2910    0.0210

% michol:
%         0         0         0         0
%    0.0205    0.1653   28.0000   10.0000
%    0.0377    0.0432   42.0000   21.0000
%    0.0592    0.0635   82.0000   21.0000
%    0.0661    0.0842  105.0000   21.0000
%    0.2047    0.1965  167.0000   21.0000
%    0.8195    0.8076  244.0000   21.0000
%    3.3675    4.1005  352.0000   21.0000
%   15.2635   24.5789  506.0000   21.0000

%no preconditioner
%   1.0e+03 *
%
%         0         0         0         0
%    0.0001    0.0001    0.0540    0.0100
%    0.0001    0.0001    0.1050    0.0210
%    0.0001    0.0002    0.2120    0.0210
%    0.0001    0.0004    0.4420    0.0210
%    0.0003    0.0011    0.9020    0.0210
%    0.0008    0.0046    1.8300    0.0210
 %   0.0034    0.0301    3.7810    0.0210
%    0.0154    0.2586    7.7280    0.0210


%		tic()
%		'lu'
%		p = symamd(AmuB);
%		lu(AmuB(p,p));
%		toc()

%{
		tic()
		nc=0;
% this does not converge
%		L = ichol(B);
		L = speye(size(B));
		[v e flag] = eigs(@(x) Afunc(B,A-mu*B,x,L),size(A,2),1,'LM',opts);
		1./e + mu
		nc
		flag
		%v'*A*v/(v'*B*v)
		%condest(A-mu*B)*condest(B)
		toc()
%}
		end
	plot(T_)

function x = Afunc(AmuB,B,x,L)
	tol = 1e-4;
	maxit = length(x);
	y = B*x;
%[x flag] = minres(AmuB,y,tol,maxit);
	x0 = []; %y./diag(AmuB);
	[x flag relres iter] = pcg(AmuB,y,tol,maxit,L,L',x0);
	if (0 ~= flag)
		'no convergence'
	else
		nc=nc+1;
		mc = mc+iter;
	end
end

end
