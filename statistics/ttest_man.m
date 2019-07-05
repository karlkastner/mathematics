% Wed 12 Dec 15:07:10 CET 2018
%
%% two-sample t-test
%% unequal sample size
%% equal variance
function [reject,p,ci,stat] = ttest_man(x1,x2)
	alpha = 0.05;
	% number of samples
	n1 = length(x1);
	n2 = length(x2);
	% sample means
	m1 = mean(x1);
	m2 = mean(x2);
	% sample variances
	sx12 = var(x1);
	sx22 = var(x2);
	% degrees of freedom
	dof = n1+n2-2;
	% pool the variance
	sp2 = ((n1-1)*sx12 + (n2-1)*sx22)/dof;
	% squared t-statistic
	t2 = (n1*n2)/(n1+n2)*(m1 - m2).^2*1/sp2;

	t  = sqrt(t2);
	% te      = tinv(1-0.5*(1-alpha),dof);
	% reject = (t>te); 
	% tcdf(t) > 1-0.5*(1-alpha) = 1/2+1/2 alpha
	p = 2*tcdf(t,dof)-1;
	reject = p>alpha;
%	p = tcdf(t,dof);

	ci = [];
	stat.tstat = t;
	stat.df    = dof;
	stat.se    = sp2;
end % ttest_man

