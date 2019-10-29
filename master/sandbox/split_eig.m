function [A B] = split_eig(n)
%{
spectral
fft -> eig -> ifft
symmetric matrix
real eigenvalues
higher accuracy
faster (less thorughput, better separated)
only positive eigenvalues
%}

clf
%n=17
%n=4097
%n=1025
%n=512
A = full(spdiags(ones(n,1)*[1 -2 1], -1:1,n,n));
B = abs(diag((0:n-1)-n/2+0.5));
As = A(1:n/2+0.5,1:n/2+0.5);
Bs=B(1:n/2+0.5,1:n/2+0.5);
tic
E1 = sort(eig(A+B));
%E1 = sort(eigs(sparse(A+B),n-1,'LM'));
t1=toc
As(end,end-1) = 0;
tic
E2 = sort(eig(As+Bs));
%E2 = sort(eigs(sparse(As+Bs),n/2+0.5-1),'LM');
t2=toc
As(end,end-1) = 2;
tic
E3 = sort(eig(As+Bs));
%E3 = sort(eigs(sparse(As+Bs),n/2+0.5-1),'LM');
t3=toc

plot(E1,'*','Markersize',10)
hold on
plot(0:2:2*length(E2)-1, E2,'or','Markersize',10)
plot(1:2:2*length(E2), E3,'og','Markersize',10)
grid on

[[NaN; E1] sort([E2; E3]) [E2; E3]];
[t1 t2+t3 t2 t3]


E23 = sort([E2; E3]);
if(length(E23) > length(E1))
	E1 = [0; E1];
end

idx=find(abs(E1 - E23) > 1e-7);
abs(E1 - E23)
[E1(idx) E23(idx)]
[E1 E23]

A
B
end % function split_eig



