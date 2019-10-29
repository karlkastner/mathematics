% Thu Sep 22 22:23:19 MSD 2011
% Karl Kästner, Berlin
%
% function [A h] = hydrogen_boxed(n, dimension)
%
% set up finite difference matrix representing a boxed hydrogen atom
%
% n : number of grid points per axis n >= 2
% dimensions : number of dimensions d \elem { 1, 2, 3 }
% A : FDM matrix
% h : grid-cell width
% x0 : location of hydrogen atom
function [A X L V h] = hydrogen_boxed(n, x_0, L0, dimension, p_fdm, p_grid, singularity_fix, mu)
	% domain size
	domain = [-L0 L0];

	% discretise domain
	% width of grid cell
	h = (domain(2)-domain(1))/(n+1);

	% a single grid-axis
	% TODO - different X for different domains
	[X_bc scale] = xgrid(n+2,p_grid,x_0,L0);
	X_bc = X_bc';
	X = scale*X_bc(2:end-1);
	%X = domain(1)+h:h:domain(2)-h;
	
	% one dimensional laplacian
%	L = laplacian(X_bc,p_fdm);
%	L = scale^-2*L;
	L = laplacian_simple(n, p_fdm)/(2*L0)^2;

	% one dimensional Coulomb potential
	% V = spdiags(1./(X - x_0)'.^2, 0, n, n);
	V = potential(X, x_0, singularity_fix, mu);
	%V=diag(1./diag(V.*V));
	V = diag(sparse(1)./diag(V)).^2;

	% step up to higher dimensions
	I = speye(n);
	%  m=3; L=full(spdiags(ones(m,1)*[1 -2 1],-1:1,m,m)); n=3; LL=0; for idx=0:n-1; LL = LL + kron(eye(m^idx),kron(L,eye(m^(n-1-idx)))), end
	switch (dimension)
		case { 1 }
			% nothing to do
		case { 2 }
			L = kron(I,L) + kron(L,I);
			% dx^2 + dy^2
			V = kron(I,V) + kron(V,I);
		case { 3 }
			L = kron(I, kron(I, L)) + kron(I, kron(L, I)) + kron(L, kron(I, I));
			% dx^2 + dy^2 + dz^2
			V = kron(I, kron(I, V)) + kron(I, kron(V, I)) + kron(V, kron(I, I));
		otherwise
			'dimension must be 1, 2 or 3'
			return
	end % switch

	% radius = sqrt sum d^2
	V = diag(sparse(1)./sqrt(sparse(diag(V))));

	% construct PDE, e.g. combine matrices to represent the non-dimensional Schrödinger equation
	% Schrödinger equation
	% ( \Laplace + V(r) ) \psi = 0
	A = 0.5*L + V;
	L = 0.5*L;
end % hydrogen_boxed

