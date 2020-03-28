% Tue Feb 28 22:55:10 MSK 2012
% Karl Kästner

function test_fem(d,mode)
opengl neverselect

k=10;
switch (d)
case {1}
	fem = @fem_1d;
	int = @int_1d_cp
	N=2.^(2:16) + 1
	eig_mode = 'SM'
	switch (mode)
		case {1}
			% discrete Laplacian
			L0 = 1
			c = 0
			f = 'f_const'
		case {2}
			% harmonic oscillator
			L0 = 10
			c = 0
			f = 'f_harmonic'
		case {3}
			% hydrogen atom
			L0 = 80
			c = L0
			%s = -100; k=1; eig_mode = 'SM'; %SA
			s = -1; eig_mode = 'SM'
			f = @f_coulomb
	end % switch mode
case {2}
	fem = @fem_2d;
	eig_mode = 'SM'
	% number of grid points
	N=2.^(3:9) + 1;
	% integration rule
	%int = {'int_2d_cp' } %, 'int_2d_smp', 'int_2d_gauss'} %, 'int_2d_corner'};
	int = @int_2d_cp;
	switch (mode)
		case {1}
			L0 = 1;
			c = 0; % zero potential
			f = @f_const;
			k_ = ceil(sqrt(k));
			E_true = pi^2*(1:k)'.^2*ones(1,k);
			E_true = E_true + E_true';
			E_true = sort(E_true(:));
			E_true = E_true(1:k);
			%E_true = pi^2*[2 5 5 8 10 10 13 13 17 17 18 20 20]';
			s = 2*pi^2;
		case{2}
			L0 = 10;
			c = L0;
			f = @f_harmonic_oscillator;
			E_true = 2*[1; 2*ones(2,1); 3*ones(3,1); 4*ones(4,1); 5*ones(5,1)];
			s = 0;
		case{3}
			L0 = 10;
			c = L0;
			f = @f_coulomb;
			s = -4.01; %2.1;
	end % switch mode
case {3}
	fem = @fem_3d;
	% number of grid points
	N=2.^(2:5) + 1;
	% integration rule
	int = {'int_3d_cp'} %, 'int_3d_smp', 'int_3d_gauss', 'int_3d_corner'};
	int = @int_3d_cp
	eig_mode = 'SM'
	switch (mode)
	case {1}
		L0 = 1;
		c = 0; % zero potential
		f = @f_const;
		k_ = ceil(sqrt(k));
		E_true = pi^2*(1:k)'.^2*ones(1,k);
		E_true = E_true + E_true';
		E_true = sort(E_true(:));
		E_true = E_true(1:k);
		E_true = pi^2*[3 6 6 6 9 9 9 12]';
	case{2}
		L0 = 10;
		c = L0;
		f = @f_harmonic_oscillator;
		E_true = [1 5 5 5 7 7 7]';
	case{3}
		L0 = 20;
		c = L0;
		f = @f_coulomb;
		E = -[    1
			1/4*ones(4,1)
			1/9*ones(9,1)
		     ];
		s = -1; % -1;
	end % switch % mode
end % switch d

E = zeros(k,length(N));

for idx=1:length(N);
	N	
	n=N(idx);
	tic()
	switch (d)
	case {2}
		% generate mesh
		[P T] = mesh_2d_uniform(n,L0);
		% set up discretisation matrices
		[A B] = feval(fem,P,T,int, @(q) feval(f,q, c));
	end
	TA(idx) = toc()
	tic()

	E(1:min(k,n^d),idx) = sort(eigs(A-s*B,B,min(k,n^d),eig_mode))+s;
	fdx=find(abs(E)>1000);
	E(fdx)=0;
	E
	
	% todo - eigenvector and error estimate
	[v e] = eigs(A-s*B,B,1,eig_mode);
	TE(idx) = toc()
	e = e+s;
	
	tic
	BC = [];
	[M err dV C] = fem_2d_mark(P, T, BC, v);
	Tr(1,idx) = toc()

	tic
	clf
	subplot(2,2,1)
	display_2d(P,T,dV(:,2)/max(dV(:,2)));
	axis equal; axis tight
	subplot(2,2,2)
	display_2d(P,T,dV(:,1)/max(dV(:,1)));
	axis equal; axis tight
	subplot(2,2,3)
	display_2d(P,T,err/max(err));
	axis equal; axis tight
	colorbar
	drawnow
	Td(1,idx) = toc()

	%E(:,idx)*(n-1)/pi^2
end % for idx

%E_true = (1:k).^2'*pi.^2;
E_true = E(:,end);
%[E_true E]
Err = E - E_true*ones(1,length(N));
nErr = sqrt(sum(Err.^2))
loglog(N,nErr,'.-')

end % function test_fem


