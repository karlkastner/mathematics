L0 = 1;
e_true = 1/L0^2*pi^2;
N = [1:100]; %1:20;
for idx=1:length(N)
	n = N(idx);
	h = L0/(n+1);
	A = 1/h^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);          
	e = eigs(A,[],1,'SM');
	err(idx,1) = abs(e+e_true);
	A_ = A - 1/12*h^2*A^2;
	e = eigs(A_,[],1,'SM');
	err(idx,2) = abs(e+e_true);
	A_ = A_ + 1/90*h^4*A^3;
	e = eigs(A_,[],1,'SM');
	err(idx,3) = abs(e+e_true);
	A_ = A_ - 1/560*h^6*A^4;
	e = eigs(A_,[],1,'SM');
	err(idx,4) = abs(e+e_true);
	A_ = A_ + 1/3150*h^8*A^5;
	e = eigs(A_,[],1,'SM');
	err(idx,5) = abs(e+e_true);
end
loglog(N,err)

	spy(A_)
