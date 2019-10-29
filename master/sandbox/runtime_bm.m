 N=round(logspace(1,4,10)); for ndx=1:length(N); n=N(ndx),
A=spdiags(ones(n,1)*[1 -2 1],-1:1,n,n)+sparse(diag(rand(n,1))); tic; m=0;
while(1); eig(A); m=m+1; t=toc; if (t>1); break; end, end, T(ndx)=toc/m; end;
loglog(N,T); grid on, set(gca,'minorgrid','none'), T

