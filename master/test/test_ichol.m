N=60:10:100;
T=[];
for idx=1:length(N)
	n=N(idx); %25;
	A=(n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n) + 0/(n+1).^2*spdiags((1:n).^2',0,n,n);
	I=speye(n);
	A=kron(kron(A,I),I)+kron(kron(I,A),I)+kron(I,kron(I,A));
	b = rand(n.^3,1);
	
	tic
	pcg(-A,-b,[],n^3);
	T(idx,1) = toc;
	
	tic
%	PC=sqrt(inv(-diag(diag(A))));
	%PC=cholinc(-A,'0');
%	PC=ichol(-A);
	PC = tril(-A) + 0.4*diag(diag(A));
	T(idx,2) = toc;
	tic
	pcg(-A,-b,[],n^2,PC,PC');
	T(idx,3) = toc
	T(idx,4) = T(idx,2)+T(idx,3)
end	

plot(N,T)
legend('CG no P','P','CG with P', 'P and CG')

