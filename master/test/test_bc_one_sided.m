[dF dF_num dF_den] = derive_bc_one_sided();
%[1.0000   -4.0000   -6.0000   20.0000  -11.0000]/12
N = 2.^(4:16);
A_ = dF_num(2,2:5)/dF_den(2);
pause
for idx=1:length(N)
	n = N(idx);
	A = spdiags(ones(n,1)*[1 -2 1], -1:1, n, n);
	E1(:,idx) = sort(eigs(A,[],10,'SM'))
	A = spdiags(ones(n,1)*[-1 16 -30 16 -1]/12, -2:2, n, n);
	E2(:,idx) = sort(eigs(A,[],10,'SM')) % symmetric, but just 2nd order
	%A(1,1:4) = A_;			% unsymmetric, 4th order, late convergence
	%A(end,end-3:end) = fliplr(A_);
	A(1,1) = A(1,1) +1/12;		% symmetric + 4th order - but only without potential
	A(end,end) = A(end,end) +1/12;
	A = A*(n+1)^2;
%	A = A + diag(sparse(1./((1:n)/(n+1))));
	E3(:,idx) = sort(eigs(A,[],10,'SM'))
%	full(A)
%	pause
end

D1 = E1 - E1(:,end)*ones(1,length(N))
D2 = E2 - E2(:,end)*ones(1,length(N))
D3 = E3 - E3(:,end)*ones(1,length(N))

D1 = sum(D1.*D1)
D2 = sum(D2.*D2)
D3 = sum(D3.*D3)
loglog(N,[D1; D2; D3]')
legend('2nd','4th no bc','4th with bc')

