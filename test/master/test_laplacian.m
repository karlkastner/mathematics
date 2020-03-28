
n=1024;
%X=(0:n+1)'/(n+1);
X=((0:n+1)'/(n+1) - 0*0.5);
X = pi*X;
h=1;
%h = 1/(n-2);
h=X(2)-X(1);
%h=1; %sqrt(34.5742)/n;
L2_test = laplacian(X,2);
L2 = 1/h^2*(spdiags(ones(n,1)*[1 -2 1], -1:1, n,n));

norm(full(L2-L2_test))

	L4_test = laplacian(X,4);

	kernel = [-1 16 -30 16 -1]/12;
        ip     = [ 5  -10   10   -5   1 ];
	L4=spdiags(ones(n,1)*kernel/h^2, -2:2, n,n);
	% boundary condition, interpolate f(-h)
	% interpolation kernel
	ik = kernel(1)/h^2*ip(2:end);
	L4(1,1:4) = L4(1,1:4) + ik;
	L4(end,end-3:end) = L4(end,end-3:end) + fliplr(ik);

	%ik =  1/12*1/h^2*[-1];
	%full(L4(1,1:3))
	%L4(1,1) = L4(1,1) - ik;
	%full(L4(1,1:3))
	%L4(1,1:3) = 1/h^2[-30 16 -1]/12;

	% linear interpolation
	%L4(1,1:3) = 1/h^2*[-29 16 -1]/12;
norm(full(L4-L4_test))

L4 = L4_test;
X_ = X(2:end-1);
b = sin(X_);
C = [X_ b L2*b L4*b L2*b+b L4*b+b L2*b-L4*b];
%[ C(:,3)./C(:,2) [diff(X_); 0]]
sqrt(sum(C.*C)) 
%[X']
%[-1+h:h:1-h+1e-7]
%h
%full(L2)
return


% laplacian test case
n=100; h=pi/(n+1);
A=(1/h^2)*full(spdiags(ones(n,1)*[1 -2 1],-1:1,n,n));
x=(1:n)'*h;
y=sin(x);
%y = x.*x.*x;
subplot(2,1,1)
plot(x,[y A*y A*A*y])
xlim([0 pi])

subplot(2,1,2)
L = spdiags(ones(n,1)*(1.0/(12*h*h))*[-1 16 -30 16 -1],-2:2,n,n); % check for correctness
%L = spdiags(ones(n,1)*(1.0/(3*h*h))*[1 2 -8 2 1],-2:2,n,n); % check for correctness
plot(x,[y L*y L*L*y])
xlim([0 pi])
%ylim([-2 2])
[y L*y L*L*y]


