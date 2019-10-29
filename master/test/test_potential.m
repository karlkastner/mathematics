
n=100; x0=1; X=2*(1:n)/n; [XX YY] = meshgrid(X,X); V=(1./sqrt(((XX-x0).^2 + (YY-x0).^2))); imagesc(V); axis square   

singularity_fix 
V = potential(X, x0, singularity_fix, mu);


