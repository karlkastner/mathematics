%plot([-1 1 1], [0 0 0],'.-')
%text( -1, 0, '1')
%text(0, 0, '-2')
%text( 1, 0, '1%')
clf
n=5;
I = eye(n); A = spdiags(ones(n,1)*[1 -2 1], -1:1, n,n);
subplot(2,3,1);
spy(A,'k');
subplot(2,3,2);
spy(kron(A,I) + kron(I,A),'k')  
subplot(2,3,3);
spy(kron(kron(A,I),I) + kron(kron(I,A),I) + kron(kron(I,I),A),'k')  
print -deps fdm_structure.eps

