% Thu  2 Jul 17:34:51 +08 2020
%p=0.5;
%=2;
% l=3;
%	mu =  a=mean(x.^p), 
%	mu = l.^p*gamma(n+p)/gamma(n)
abc = perms([2,3,5]);
mu  = [];
n   = 1e6;
for idx=1:size(abc,1)
	a = abc(idx,1);
	b = abc(idx,2);
	c = abc(idx,3);

	x=gamrnd(a,b,n,1);
	%x=gamrnd(n,l,1e6,1);
	mu(idx,1) = mean(x.^c);
	mu(idx,2) = mean_generalized_gampdf(a,b,c);
end
disp(mu)
mu./mu(:,1)
