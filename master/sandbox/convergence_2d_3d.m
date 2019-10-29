% 2012 Jan 23 14:01
% Karl Kästner, Berlin

% inifinite space, nucleus in the corner (symmetrie)

function convergence_2D_3D()

L0 = 1000;
N = 2.^(2:10);

for ndx=1:length(N);
	n = N(ndx);
	A = setup_3D(L0,n);
	I2 = speye(size(A));
	s = -0.5; %-0.0808;
	E(:,ndx) = sort(eigs(A-s*I2,20,'SM')) + s

	c = derive_richardson(ndx)';
	F(:,ndx) = E(:,1:ndx)*c;
	1./F
end

end % convergence_2D_3D

function A = setup_3D(L0,n)
	X1 = sparse(L0*(1:n)'/(n+1));
	I1 = speye(n);
	L1 = (n+1)^2/L0^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	L  = kron(I1,kron(I1, L1)) + kron(I1,kron(L1, I1)) + kron(I1,kron(I1,L1));
	Rs =   kron(I1,kron(I1,diag(X1.^2))) ...
	     + kron(I1,kron(diag(X1.^2),I1)) ...
             + kron(I1,kron(I1,diag(X1.^2)));
	V  = diag(sparse(1)./diag(sqrt(Rs)));
	A = -0.5*L - V;
end

function A = setup_2D(L0,n)
	X1 = sparse(L0*(1:n)'/(n+1));
	I1 = speye(n);
	L1 = (n+1)^2/L0^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	L  = kron(I1, L1) + kron(L1, I1);
	Rs = kron(I1,diag(X1.^2)) + kron(diag(X1.^2),I1);
	V  = diag(sparse(1)./diag(sqrt(Rs)));
	A  = -0.5*L - V;
end

