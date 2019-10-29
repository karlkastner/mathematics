% Sun Feb 26 15:45:25 MSK 2012

function E = test_fem_3d(mode)

% number of grid points
N=2.^(2:5) + 1;
% number of eigenvalues
k = 10;
% numer of eigenvaluess to compare
ke = 10;

int = {'int_3d_cp'} %, 'int_3d_smp', 'int_3d_gauss', 'int_3d_corner'};
colour = {'b','r','k','g'};
clf

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
	s = -100;
case{2}
	L0 = 10;
	c = L0;
	f = @f_harmonic_oscillator;
	E_true = [1 5 5 5 7 7 7]';
	s = -2;
case{3}
	L0 = 20;
	c = L0;
	f = @f_coulomb;
	E = -[    1
		1/4*ones(4,1)
		1/9*ones(9,1)
	     ];
	s = -1; %0.5;
end % switch

for idx=1:length(int)

E=zeros(k,length(N));
for ndx=1:length(N)
	N
	n=N(ndx)
	% matrix setup
	tic
	[A B] = fem_3d(n, L0, int{idx}, @(q) f(q,c));
	TA(ndx) = toc()
	% solution of the eigenvalue problem
	tic
	E(1:min(k,n^3),ndx)=sort(eigs(A-s*B,B,min(k,n^3),'SM'))+s;
	E(1:ke,:)
	[bcdx] = find(E>1e3);
	E(bcdx)=0;
	TE(ndx)=toc()
end % ndx

% verify the result

E_true = E(:,end);
[E_true(1:ke,:) E(1:ke,:)]
Err = E(1:ke,:) - E_true(1:ke,:)*ones(1,length(N));
nErr = sqrt(sum(Err.^2))
loglog(N,nErr, [colour{idx} '.-'],'LineWidth',2,'Markersize',20);
hold on
grid on

end % for idx

legend(int{:});

end % test_fem_3d

