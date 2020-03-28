% Tue Jun  5 00:16:06 MSK 2012
% Karl Kästner, Berlin

function test_preconditioning

tol = sqrt(eps);
n = 200;
maxit = n^2;
% get the 2D poison problem
	I = speye(n); 	A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n); 	A = kron(I,A) + kron(A,I);
% get an input vector
	b = rand(n^2,1);

% invert signs
A = -A;
b = -b;

% get a preconditioner
L = ichol(A);

% solve without precondtioner
tic()
minres(A,b,tol,maxit,[],[],zeros(n^2,1));
toc()

% solve with precondtioner
tic()
minres(A, b, tol, maxit, L', L, zeros(n^2,1));
toc()

% solve with precondtioner function
tic()
minres(A, b, tol, maxit, @(x) lfunc(L',x), @(x) lfunc(L,x), zeros(n^2,1));
toc()

tic()
eigs(A,[],1,'SM')
toc()

tic()
eigs(@(x) afunc(A,x,L),size(A,1),1,'SM')
toc()


function x = lfunc(L,x)
	x = L \ x;
end

function x = afunc(A,x,L)
	%x = A \ x;
	%x = minres(A,x,tol,maxit);
	%x = minres(A,x,tol,maxit,L',L);
	x = pcg(A,x,tol,maxit,L',L,zeros(size(x)));
end

end

