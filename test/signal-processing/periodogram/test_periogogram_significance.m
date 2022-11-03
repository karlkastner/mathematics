n=1e6; x=randn(n,1); S=periodogram(x,1); p=(1:n)'/(n+1); plot([sort(S),2/n*sort(1/2*chi2rnd(2,length(S),1)),1/n*chi2inv(p,2)])
