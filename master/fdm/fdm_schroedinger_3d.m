% Wed May 23 23:12:19 MSK 2012
% Karl Kästner, Berlin

% Sets up the variable grid FDM discretisation matrix
% of the Schrödigner equation in 3D
function [A B D] = fdm_schroedinger_3d(X)
	X1 = X{1};
	X2 = X{2};
	X3 = X{3};
	n1 = size(X1,1)-2;
	n2 = size(X2,1)-2;
	n3 = size(X3,1)-2;
	I1 = speye(n1);
	I2 = speye(n2);
	I3 = speye(n3);
	% dummy mass matrix
	B = speye(n1*n2*n3);
	% Laplacian setup
	[L1 D1] = fdm_d_vargrid(X1,2,3);
	[L2 D2] = fdm_d_vargrid(X2,2,3);
	[L3 D3] = fdm_d_vargrid(X3,2,3);
	% combine matrices with similarity transform to preserve symmetrie
	L = kron(kron(I1,inv(D2)),inv(D3))*kron(kron(L1,I2),I3)*kron(kron(I1,D2),D3) ...
	  + kron(kron(inv(D1),I2),inv(D3))*kron(kron(I1,L2),I3)*kron(kron(D1,I2),D3) ...
	  + kron(kron(inv(D1),inv(D2)),I3)*kron(kron(I1,I2),L3)*kron(kron(D1,D2),I3);
	% matrix to undo symmetry transform
	D = kron(kron(D1,D2),D3);
	% Coulomb potential
	% strip boundary points
	R1 = abs(X1(2:end-1));
	R2 = abs(X2(2:end-1));
	R3 = abs(X3(2:end-1));
	R  = sqrt( kron(kron(R1.^2, ones(n2,1)),ones(n3,1)) ...
		 + kron(kron(ones(n1,1), R2.^2),ones(n3,1)) ...
		 + kron(kron(ones(n1,1), ones(n2,1)), R3.^2) );
	V  = diag(sparse(1./R));
	% Schrödinger equation
	A = -(0.5*L + V);
end % fdm_schroedinger_3d

