% Tue Jun 19 18:13:28 MSK 2012
% Karl Kästner, Berlin

function example_2d(mode)
	javaaddpath('.');
	javaaddpath('/usr/share/java/jama.jar');
	switch(mode)
	    case {1}
		% discrete Laplacian
		k=1;
		opt.shift = 1.1*pi^2;
		E_true = pi^2*(1:k)'.^2*ones(1,k); E_true = E_true + E_true'; opt.E_true = sort(E_true(:));
		pdeeig_2d(1, [], 1, [1 1], [], opt);   
	    case {2}
		% harmonic oscillator
		opt.shift = 1.1*2;
		opt.E_true = 2*[1; 2*ones(2,1); 3*ones(3,1); 4*ones(4,1); 5*ones(5,1)];
		pdeeig_2d(1, @f_harmonic_oscillator, 1, 100*[1 1], [], opt);   
	    case {3}
		% hydrogenic Schrödinger equation
		opt.shift = -2.2;
		opt.E_true = -0.5./([1 2*ones(1,3) 3*ones(1,5) 4*ones(1,7) 5*ones(1,9) 6*ones(1,11) 7*ones(1,13)]-1/2).^2; % number : 2*n - 1
		opt.reltol = 1e-4;
		opt.backward = 1;
		opt.circular = 1;
		opt.order = 4;
		k = 1;
		%L0 = 120;
		L0 = 1;
		f = Potential_2D_Coulomb;
		%opt.order = 2; L0 = 10; opt.mflag = 2;   
		pdeeig_2d(-0.5, f, 1, L0*[1 1], [], opt);   
	end % switch mode
end % example_2d()


