L0=0.5; p=3; n=1000; X=(0:n)/(n) - 0.5; X=L0*sign(X.^(p-1)).*X.^p/X(1)^p; A=laplacian_non_uniform(X'); eigs(A,[],5,'SM')/pi^2

