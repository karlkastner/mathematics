% Jul  9 23:32 (MSK)
% Fri Jul 13 21:09:17 MSK 2012
% Karl Kästner, Berlin

syms l h a
f_ = (l*4*h/6 - 2/h)*sin(a*h) + (l*h/6 + 1/h)*sin(a*2*h)
f_ = solve(f_,l)
pause

N=1:100;
for idx=1:length(N)
	n=N(idx);
	M = 1/(6*(n+1))*spdiags(ones(n,1)*[1 4 1],-1:1,n,n);
	K = (n+1)*spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);
	A = (n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n);
	if (1 == n)
		err(idx,1) = eig(-full(K),full(M)) + pi^2;
		err(idx,2) = eig(full(A)) + pi^2;
	else
		err(idx,1) = eigs(-K,M,1,'SM') + pi^2;
		err(idx,2) = eigs(A,[],1,'SM') + pi^2;
	end
	L = chol(M);
	6*full(M)*(n+1)
	full(K)/(n+1)
	full(L)
	full(inv(L))
	full(inv(L')*K*inv(L)/(n+1)^2)
	pause
end
err = abs(err);

%loglog(N,err)
plot(N,err)


