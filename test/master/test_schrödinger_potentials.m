
two regions
p=1e3; k=5; n=10000; A = -(n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n) + diag(sparse( p*(((1:n)-n/2) < 0)) ); [v e] = eigs(A, k, 'SM'); plot(v.^2)     
spike in centre
p=1e4; k=5; n=10000; A = -(n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n) + diag(sparse( p*( abs((1:n)-n/2) <= 1)) ); [v e] = eigs(A, k, 'SM'); plot(v.^2)
% harm osci
p=1e3; k=5; n=10000; A = -(n+1)^2*spdiags(ones(n,1)*[1 -2 1],-1:1,n,n); D=diag(sparse( 10*(((1:n)-n/2)/(n+1))).^2 ); [v e] = eigs(A+D, k, 'SM'); plot(v.^2)  

