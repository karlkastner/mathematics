% Mon Oct 10 23:18:09 MSK 2011
% Karl Kästner, Berlin

% n number of grid points in r, phi and z axis
% r radius
% h height
function D2 = laplacian_cylindrical(n, r, h)
	% accept only even grids
	if (0 ~= bitand(n,1))
		'error'
	end
	% stepwidth
	dphi = 2*pi/n;
	dr = r/n;
	dz = h/n;

	I = speye(n);
	R = spdiags(((0:n-1)+0.5)*dr,0,n,n);
	dR  = 1/dr*spdiags(ones(n,1)*[-1 0 1],-1:1,n,n);
	BCdR = sparse(1,1,(1/dr)/R(1),n,n);
	ddR = 1/dr^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	S = spdiags(ones(n/2,1),-n/2,n,n);
	BCddR = sparse(1,1,1/dr^2,n,n);

	% Phi-axis
	ddPhi = 1/dphi^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	% periodic BC
	ddPhi(1,end) = 1/dphi^2;
	ddPhi(end,1) = 1/dphi^2;

	% Z-axis with dirchlet BC - simple laplacian
	ddZ = 1/dz^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	
	% 3D
	% L = 1/r*d/dr + d^2/dr^2 + 1/r^2*d^2/dphi^2 + d^2/dz^2
	L =       kron(kron((1./R)*dR,I),I) ...
		+ kron(kron(S,I),BcdR) ...
		- kron(kron(S',I),BCdR) ...
                + kron(kron(ddR,I),I) ...
		+ kron(kron(S,I),BCddR) ...
		+ kron(kron(S',I),BCddR) ...
		+ kron(kron(1./R.^2, ddPhi), I) ...	% todo, check, if correct axis scaled
		+ kron(kron(I,I),ddZ);
end % function laplacian_cylindrical

