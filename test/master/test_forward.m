%For standard laplacian without potential:
%test 1: just one sided FDM close to x0
%test 2: cut out region around x0 and extrapolate

k = 10;

kc = [1 -2  1];

N=2.^(4:8);

for idx=1:length(N);
	n = N(idx);
	h = 1/(n+1);
	x = 2*(((1:n)/(n+1)) - 0.5);
	X = diag(sparse(x));
	A = 1/h^2*spdiags(ones(n,1)*kc,-1:1,n,n);
	B = (kron([1 0; 0 0],spdiags(ones(n/2,1)*kc, -2:0, n/2,n/2)) + kron([0 0; 0 1], spdiags(ones(n/2,1)*kc,  0:2, n/2,n/2)));
	B(1,:) = 0; B(1,1:3) = [-5 4 1]; B(end,:) = 0; B(end,end-2:end) = [1 5 -5];
	B = 1/h^2*B;

	r = 2*fliplr(((1:n/2)/(n+1)))'
	rr = 2*(((1:n/2)/(n+1)))'
	A = (kron([1 0; 0 0],spdiags([(2-r) -2*ones(n/2,1) r], -1:1, n/2,n/2))) + kron([0 0; 0 1], spdiags([rr -2*ones(n/2,1) (2-rr)],  -1:1, n/2,n/2));
	A = A/h^2;
%full(A)
%pause

	I = speye(n);
	AA = kron(A,I) + kron(I,A);
	BB = kron(B,I) + kron(I,B);
	R = sqrt(kron(I,X.^2) + kron(X.^2,I));
	R_ = diag(sparse(1-diag(R)));
	CC = R*AA + R_*BB;
%	full(diag(R))
%	full(diag(R_))
%	pause
%	C = diag(sparse(R))*A + diag(sparse(1-R))*B;
%	E(:,idx) = sort(eigs(AA,k,'SM'));
%	E(:,idx) = sort(eigs(BB,k,'SM'));
%	E(:,idx) = sort(eigs(B,k,'SM'));
	E(:,idx) = sort(eigs(A,k,'SM'));
end

E
E = E - E(:,end)*ones(1,length(N))

