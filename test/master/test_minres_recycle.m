n=3;
d=2;
k = 10;
%A=rand(n); A=A+A';
%A = rand(n^2); A=A*A' + n^2*speye(n^2); %A=A+A';
%A = A; %- n^2*speye(n^2); %diag(sparse(1:n^2)).^2;
A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n); I = speye(n);
switch (d)
	case {2}
		% Poisson 2D
		A = kron(A,I) + kron(I,A);
	case {3}
		% Poisson 3D
		A = kron(kron(A,I),I) + kron(kron(I,A),I) + kron(I,kron(I,A));
end
% make spd
A=-A;
% k initial vectors
b = rand(n^d,k);

% solve repetitively the k linear systems
clf
for idx=1:k
	[x V T r beta R R2] = minres_recycle(A,b(:,idx),1e-7,[],[],4*n^d);
%	V*T - A*V
%	x - (A\b(:,idx))
%	full(A)
%	V*T*V' - A
%	[x A \ b(:,idx) ]
	norm(x - A \ b(:,idx))
	subplot(2,2,1)
	semilogy(1:length(beta),beta); hold on	
	ylim([1e-15 1e15])
	subplot(2,2,2)
	semilogy(1:length(R),R); hold on	
	subplot(2,2,3)
	semilogy(1:length(R),R2); hold on	
end % for k
%hold off
