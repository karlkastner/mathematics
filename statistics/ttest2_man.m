% Sun Oct 13 16:36:18 UTC 2013
% Karl Kästner, Berlin
%
%% two-sample t-test
%% here posix return value standard: h = 0 accepted, h = 1 failed
%% note: the matlab logic is inverse : h = 1 accepted, h = 0 failed
%% two sided univariate t-test
%
% alpha : confidence level
function [reject,alpha_e,ci,stat] = ttest2_man(x1,x2,alpha, vflag)
	if (nargin()<3 || isempty(alpha))
		alpha = 0.05;
	end
	if (nargin() < 4)
		vflag = true;
	end

	% number of samples
	n1 = length(x1);
	n2 = length(x2);

	% sample means
	x1_bar = mean(x1);
	x2_bar = mean(x2);
	% sample variances
	% note : original code missed normalization by n1-1
	s12 = 1/(n1-1)*sum((x1 - x1_bar).^2);
	s22 = 1/(n2-1)*sum((x2 - x2_bar).^2);
	% number of degrees of freedom
	if (nargin() < 4 || ~vflag)
		% equal variances - pooled variance
		df = n1+n2-2;
		% pooled sample variances
		s2 = ((n1-1)*s12 + (n2-1)*s22 )/(df);
		se = sqrt((n1+n2)/(n1*n2)*s2);
		t2 = (n1*n2)/(n1 + n2)*(x1_bar - x2_bar)*(1/s2)*(x1_bar - x2_bar);
	else
		% unequal variances - unpooled variance
		s2 = s12/n1 + s22/n2;
		df = (s12/n1 + s22/n2)^2/( (s12/n1)^2/(n1-1) + (s22/n2)^2/(n2-1));
		t2 = (x1_bar - x2_bar)*(1/s2)*(x1_bar - x2_bar);
	end

	% t-statistic (squared)
	%t2 = (n1*n2)/(n1 + n2)*(x1_bar - x2_bar)*(1/s2)*(x1_bar - x2_bar);

	% t-statistic
	t  = sqrt(t2);
	% probability, that samples are not equal
	% te = tinv(1-0.5*(1-alpha),df);
	% reject = (t > te);
	% p-value = alpha_e
	% TODO matlab defines this as 2*tcdf(-t,df)
	alpha_e = 2*tcdf(t,df) - 1; 
	reject = alpha_e > alpha;	

	ci     = [];
	stat.dif = (x1_bar - x2_bar);
	stat.t = t;
	%stat.te = te;
	stat.df = df;
	stat.sd = sqrt(s2);
	stat.se = se;
end % ttest2_man()

