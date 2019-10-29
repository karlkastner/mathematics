% Thu Nov  3 18:19:48 MSK 2011
% Karl Kästner, Berlin

% test the computation of the energy levels of an elelctron in a quantum well
function test_quantum_well()
	% number of eigenvalues to compute
	m = 16;
	% problem setup
	N = [128];
	L = [1.0];
	X0 = []; % solely electron, no proton
	% compute the eigenvalues
	E_numeric  = .. (N,L,X0);
	E_analytic = quantum_well(N, L, X0);
	% calculate the error norm
	err = sort(E_numeric) - sort(E_analytic);
end % test quantum well

% compute energy levels of an electron in an infinite potential well
% well = box with side length 2 L a0
function E = quantum_well_levels(N, L)
	a0 = 5.2917721092e-11;  % Bohr radius in meter m
	h_bar = 6.626e-34/(2*pi); % reduced Planck's constant in J*s
	me = 9.109e-31;	% electron mass in kg

	a = 2*L*a0;
	switch (N)
		case {1}
			k = (1:n)/a(1);
		case {2}
			I1 = ones(1,N(1));
			I2 = ones(1,N(2));
			k = kron(1:N(1),I2)/a(1) + kron(I1,1:N(2))/a(2);
		case {3}
			I1 = ones(1,N(1));
			I2 = ones(1,N(2));
			I3 = ones(1,N(3));
			k = kron(kron(1:N(1),I2),I3)/a(1) ...
		          + kron(kron(I1,1:N(2)),I3)/a(2) ...
		          + kron(kron(I1,I2),1:N(3))/a(3);
	end % switch
	k = sum(k);
	E = h_bar^2 * pi^2 * k^2 / (2 * me);
end % quantum_well

