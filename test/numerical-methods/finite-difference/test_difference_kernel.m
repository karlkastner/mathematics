% Wed 26 Jun 11:24:20 CEST 2019

h = 0.1;
for d=1:4
	K=difference_kernel(d,h);
	K*h^d
end
d = 2;
o = 1;
K =	difference_kernel(d,h,o)

% the setup with matrix powers does not work at the boundary
% -> se up larger matrix, mutiply with extrapolation matrix
n=7;
D1 = 1/2*full(spdiags(ones(n,1)*[-1,0,1],-1:1,n,n));
D2 = full(spdiags(ones(n,1)*[1,-2,1],-1:1,n,n)); 2(1,1:3) = [1,-2,1]; D1(1,1:3) = [-1,1,0]; D1(end,end-2:end) = [0,-1,1]; D2(end,end-2:end)=[1,-2,1]; ;
y=((-3:n-4).^3)'; D2*D1*y, D1*D2*y

