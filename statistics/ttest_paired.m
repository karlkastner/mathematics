% Wed 12 Dec 15:07:10 CET 2018
%
%% paired t-test
%% unequal sample size
%% equal variance
%% more powerfull than unpaired test, as long as correlation between x1 and x2 > 0
function [reject,p,ci,stat] = ttest_paired(x1,x2)
	alpha = 0.05;
	n = length(x1);
	if (n~=length(x2))
		error('vectors have to be of the same length');
	end
	d    = x1-x2;
	% dbar = mean(x1-x2) = m1-m2
	dbar = mean(d);
	se   = 1/sqrt(n)*std(d)
	%s1 = std(x1);
	%s2 = std(x2);
	%se = sqrt(1/n*(s1^2+s2^2));
	t2   = dbar^2/se^2;
	dof  = n-1; 

	t  = sqrt(t2)
	te = tinv(1-0.5*(1-alpha),dof)
	reject = (t>te); 
	p = tcdf(t,dof);
	ci = [];
	stat.tstat = t;
	stat.df    = dof;
	stat.se    = se;
end % ttest_man

