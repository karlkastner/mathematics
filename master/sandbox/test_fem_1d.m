% Tue Feb 28 22:55:10 MSK 2012
% Karl Kästner

function test_fem_1d()

N=2.^(2:16) + 1
k=10;

%L0 = 1
%c = 0
%f = 'f_const'

L0 = 80
c = L0
s = -1;
k=10; s=-100;

mode = 'SM'
mode = 'SA'

E = zeros(k,length(N));

for idx=1:length(N);
	N	
	n=N(idx);
	tic()
	X = L0*(0:n)'/(n-1);
	[P T] = mesh_1d_uniform(X);
	[A B] = fem_1d(P,T,'int_1d_cp', @(x) f_coulomb(x,0));
	TA(idx) = toc()
	tic()
	C = A-s*B;
	norm(C-C')

	E(1:min(k,n),idx) = sort(eigs(A-s*B,B,min(k,n),'SM'))+s;
	fdx=find(abs(E)>1000);
	E(fdx)=0;
	E

	TE(idx) = toc()
	%E(:,idx)*(n-1)/pi^2
end

E_true = (1:k).^2'*pi.^2;
E_true = E(:,end);
%[E_true E]
Err = E - E_true*ones(1,length(N));
nErr = sqrt(sum(Err.^2))
loglog(N,nErr,'.-')

end % function test_fem_1d

