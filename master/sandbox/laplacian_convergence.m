% eigenvalues of the D-laplacian exact
% eigenvectors converge like h^2
N=2.^(0:18); for idx=1:length(N); n=N(idx); A=(n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n); [v e]=eigs(A,[],1,'SM'); err(idx,1)=abs(e/pi^2+1); v2=sin(pi*(1:n)'/(n+1)); v2=v2/norm(v2); err(idx,2)=norm(v-v2)/sqrt(n); end; loglog(N,err)

