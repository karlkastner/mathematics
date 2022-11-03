k=2; s=3; e=exprnd(1,1e6,k); e=mean(e,2); mean(e), std(e), kurtosis(e), g=1/k*gamrnd(k/s^2,s^2,1e6,1); mean(g),std(g),kurtosis(g)



