% 2015-06-22 17:45:29.219215105 +0200
% Karl Kastner, Berlin

% population size
n = 1e2;
% repetitions
m = 1e4;
% correlation
p = 0.9;
N = (1:n)';
a = sqrt((1+p)/(1-p));

title_C = { 'uncorrelated', ...
	    'uncorrelated finite', ...
	    'correlated infinite', ...
	    'correlated finite' };
	

namedfigure(1,'Correlation standard error scale factor test')
clf();
for kdx=4

e1 = randn(n,m);
if (kdx == 3 || kdx == 4)
% generate a correlated series
for idx=2:n
	e1(idx,:) = a*(1-p)*e1(idx,:) + p*e1(idx-1,:);
end
end
s = std(e1,[],2);
if (kdx == 2 || kdx == 4)
	% finite population size
	e1 = bsxfun(@minus,e1,mean(e1));
% normalise
%e1 = bsxfun(@times,e1,1./std(e1,[],2));
end
s0       = std(e1,[],2);
sd_range = mean(cumstd(e1),2);

% mean of i-successive correlated samples
mu = cummean(e1(:,:));

% standard deviation of the mean
se = std(mu,[],2);

% number of samples
n = length(se);
N = (1:n)';

%fc = f_correlation(p,N);
%f0 = sqrt((1+p)/(1-p))*x.^0;
%f1 = sqrt((1-p.^(x+1)+p)./(1-p)); %.*(1./sqrt(x))
%f1 = sqrt((1-p.^(x+1)+p)./(1-p+p.^(x+1))); %.*(1./sqrt(x))

scale = [];
scale(:,1) = se;
for idx=1:n
	switch (kdx)
		case {1} % uncorrelated infinite
			p_ = 0;
			n_ = 1e7;
			%scale(idx,2) = sqrt(ar1_var_factor(0,1e7,idx));
		case {2} % uncorrelated finite
			p_ = 0;
			n_ = n;
			%scale(idx,2) = sqrt(ar1_var_factor(0,n,idx));
		case {3} % correlated infinite
			p_ = p;
			n_ = 1e7;
			%scale(idx,2) = sqrt(ar1_var_factor(p,1e7,idx));
		case {4} % correlated finite
			p_ = p;
			n_ = n;
	end
scale(idx,2) = sqrt(ar1_var_factor(p_,n_,idx));
sigma_ = 1;
s0(idx,2) = sqrt(ar1_mse_mu_single_sample(sigma_,p_,idx,n_));
sd_range(idx,2) = ar1_var_range2(sigma_,p_,n,idx);

end

scale = bsxfun(@times,scale,sqrt(N));
%[sc ./ (1./sqrt(x)) fc f0 f1 f2])



% plot standard error for each sample (e.g. along cross section)
subplot(4,4,1+(kdx-1)*4);
plot(s0);
title('standard errors of individual samples w/r to pop mean')

% plot standard error of the mean
subplot(4,4,2+(kdx-1)*4);
plot([1./sqrt(N) sqrt((n-N)./(N.*(n-N))) se] );
title('standard error of sample mean w/r to pop mean');

%subplot(4,4,3+(kdx-1)*4);
%plot([se 1./sqrt(N)]);
%legend('se of correlated series','se of uncorrelated series');

subplot(4,4,3+(kdx-1)*4);
plot(N,scale(:,1),'linewidth',2);
hold on
plot(N,scale(:,2),'linewidth',1);
legend('simulated', 'analytic');
title('Scale factor of the standard error')
%true','complex','simple 0','simple 1')

subplot(4,4,4+(kdx-1)*4);
plot(N,sd_range)
hold on

rho(kdx) = corr(e1(end,:)',e1(end-1,:)')
err(kdx) = rmse(scale(:,1)-scale(:,2))
end % kdx

