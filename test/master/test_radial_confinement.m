n = 10240 %1
L0 = 10^3; % 10^4

%n = 102401;
%L0 = 10;
%n = 100; L0 = 10;

h = L0/(n+1);
% set up laplacian operator
D2 = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
% set up potential
X = h*(1:n);
% 
V = diag(sparse(1./X));
% compute smallest and largest eigenpair and estimate the error
A = -0.5*D2 - V;
I = speye(size(A));
s = -0.5;
%eigs(A-s*I,[],30,'SM')+s
[v e] = eigs(A,[],1,'LM');
%[v e] = eigs(A-s*I,[],1,'LM');
%e=e+s
vd4 = D2*(D2*v);
nv = sqrt(vd4'*vd4)
err = 1/12*h^2*nv

% 28 states
% analytic results for radial confinement

