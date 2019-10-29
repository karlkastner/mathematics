n = 400;
A = delsq(numgrid('L', n));
b = rand(size(A,1),1);

% A = full(A);
%b = sum(A, 2);

p = symamd(A);
tic();
[L U D P] = lu(A(p,p));
toc()

tic();
[L D P] = ldl(A(p,p));
toc()

tic();
L = chol(A(p,p));
toc()

%LSopts.POSDEF = false;
%LSopts.SYM = false;
%tic; linsolve(A, b, LSopts); toc;
%LSopts.SYM = true;
%tic; linsolve(A, b, LSopts); toc;
%LSopts.POSDEF = true;
%tic; linsolve(A, b, LSopts); toc

