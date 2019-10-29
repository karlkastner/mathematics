% Fri Feb 24 01:33:56 MSK 2012
% Karl Kästner, Berlin

function E = test_fem_2d(mode)

% number of grid points
N=2.^(3:9) + 1;
% number of eigenvalues
k = 10;
% numer of eigenvaluess to compare
ke=10;

L0 = 20;
c=L0;
f = @f_coulomb;
s = -2;

int = {'int_2d_cp' } %, 'int_2d_smp', 'int_2d_gauss'} %, 'int_2d_corner'};
colour = {'b','r','k','g'};
clf

%mode = 1

for idx=1:length(int)
	int{idx}

E=zeros(k,length(N));
for ndx=1:length(N)
	n=N(ndx)
	% matrix setup
	tic
	[A B] = fem_2d(n,L0, int{idx}, @(q) f(q,c));
	TA(ndx)=toc()
	% solution of the eigenvalue problem
	tic
	E(1:min(k,n),ndx)=sort(eigs(A-s*B,B,min(k,n),'SM'))+s
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

end % test_fem_2d

