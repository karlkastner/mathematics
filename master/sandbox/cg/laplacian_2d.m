function A = laplacian_2d(m, n)
	% n=4; m=n+1;	
	%A = rand(n*m); A=A+A';
	A1=(n+1)^2*spdiags(ones(n,1)*[1 -2 1], -1:1,n,n);
	A2=(m+1)^2*spdiags(ones(m,1)*[1 -2 1], -1:1,m,m);
	A = kron(speye(n),A2) + kron(A1, speye(m));
end

