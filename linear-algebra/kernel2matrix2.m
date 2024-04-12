% 2024-03-02 16:46:36.619430841 +0100
function A = kernel2matrix2d(a,n,dx)
	A1 = spdiags(ones(n(1),1)*(a(1,:)/dx(1)),-1:1,n(1),n(1));
	A1(1,end) = a(1,1)/dx(1);
	A1(end,1) = a(1,3)/dx(1);
	A2 = spdiags(ones(n(2),1)*(a(2,:)/dx(2)),-1:1,n(2),n(2));
	A2(1,end) = a(2,1)/dx(2);
	A2(end,1) = a(2,3)/dx(2);
	%A = kron(A1,speye(n(2))) + kron(speye(n(1)),A2);
	A = kron(speye(n(2)),A1) + kron(A2,speye(n(1)));
end

