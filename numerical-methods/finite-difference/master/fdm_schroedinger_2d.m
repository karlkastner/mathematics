% Mon May 21 23:00:54 MSK 2012
% Karl Kästner, Berlin

% Set up finite difference discretisation matrices
% for the Schrödinger equation in 2d
function [A B D] = fdm_schroedinger_2d(X)
	n1 = size(X{1},1)-2;
	n2 = size(X{2},1)-2;
	I1 = speye(n1);
	I2 = speye(n2);
	% mass matrix
	B = speye(n1*n2);
	% Laplacian setup
	[L1 D1] = fdm_d_vargrid(X{1},2,3);
	[L2 D2] = fdm_d_vargrid(X{2},2,3);
	% combine matrices with similarity transform to preserve symmetrie
	L = kron(I1,inv(D2))*kron(L1,I2)*kron(I1,D2) + kron(inv(D1),I2)*kron(I1,L2)*kron(D1,I2);
	D = kron(D1,D2);
%	D = kron(D1,I2)*kron(I1,D2);
	% Coulomb potential
	% strip boundary points
	X1 = X{1};
	X2 = X{2};
	R1 = abs(X1(2:end-1));
	R2 = abs(X2(2:end-1));
	R  = sqrt( kron(R1.^2, ones(n2,1)) + kron(ones(n1,1), R2.^2) );
	V  = diag(sparse(1./R));
	% Schrödinger equation
	A = -(0.5*L + V);
end % fdm_hydrogen_2d

