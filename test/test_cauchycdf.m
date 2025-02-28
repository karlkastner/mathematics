mu = 2;
gamma = 3;

x = [-1,1];

c = cauchycdf(x,mu,gamma)
[mu_,gamma_] = cauchycdf_fit(x,c)
