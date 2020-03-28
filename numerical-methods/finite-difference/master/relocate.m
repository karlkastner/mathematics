n=10;
A = [linspace(1,10,n), linspace(10,11,n) linspace(12,21,n)]'
n=length(A)+1;
	M=full(spdiags( [ [0.5*ones(n-3,1); 0; 0] [1; zeros(n-3,1); 1] [0; 0; 0.5*ones(n-3,1)]], -1:1, n-1,n-1))
%	M=full(spdiags( [ [0.5*ones(n-3,1); 1] [1; 0.5*ones(n-3,1)]], -1:0, n,n-1))
	size(M)
	size(A)
B(1,:) = A'
for idx=2:100
	A = M*A;
	B(idx,:) = A';
	plot(B(idx,:),ones(size(B(idx,:))),'.');
	pause(0.1);
	drawnow();
	[max(abs(diff(diff(A)))) min(diff(A)) prod([max(abs(diff(diff(A)))) min(diff(A))])]
end
