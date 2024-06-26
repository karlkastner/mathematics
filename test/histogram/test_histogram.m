% 2015-06-02 10:47:42.707088103 +0200
% Karl Kastner, Berlin

% TODO quantile implementation are not yet bias corrected
% TODO other distributions

method = {     'mixture',
	       'constant',
	       'midpoint',
	       'trapezoidal',
%	       'twopoint'
%	       'gaussian',
%	       'simpson'
	};
%method = {'mixture','constant','midpoint','trapezoidal','simpson'};

n  = 1e3;
%dh = cvec(2.^(-4:2));
dh = cvec([1/8,1/4,1/2,1,2]);
%dh = 0.5;
%dh = 0.1;
%L = 12;
L = 10;
clear bias sd bias_
for jdx=1:length(method)
for idx=1:length(dh)
	edges      = -L:dh(idx):L;
	H = Histogram([],edges);
%	H.SHEPPARD = 0;
%	H.METHOD   = method{jdx};
	H.method = method{jdx};
	
%	centres = 0.5*(edges(1:end-1)+edges(2:end));
	% draw random distribution parameter
	mu0 = 4*(rand(n,1)-0.5);
	s   = 2;
	% construct histograms
	H.setup(@(x) normcdf(1/s*bsxfun(@minus,x,mu0)));
	% estimate
	mu  = H.mean();
	s2  = H.var();
%	mean(s2)
%	var(s2)
	sk  = H.skewness();
	ku  = H.kurtosis();
	q50 = H.median();
	q16 = H.quantile(normcdf(-1));
	q05 = H.quantile(normcdf(-2));
	q84 = H.quantile(normcdf(+1));
	% bias
	bias(idx,1) = mean(mu-mu0);
	bias(idx,2) = mean(s2-s^2);
	bias(idx,3) = mean(sk-0);
	bias(idx,4) = mean(ku-3);
	bias(idx,5) = mean(q50-mu0);
	bias(idx,6) = mean(q16 - (-s+mu0));
%	bias(idx,7) = mean(q84 - (+s+mu0));
%	bias(idx,7) = mean(q05 - (-2*s+mu0));
	bias(idx,7) = mean((q84-q16)/2-s);
	bias(idx,8) = mean(q16-2*q50+q84);
	% standard deviation (not rmse!)
	sd(idx,1) = std(mu-mu0);
	sd(idx,2) = std(s2-s^2);
	sd(idx,3) = std(sk-0);
	sd(idx,4) = std(ku-3);
	sd(idx,5) = std(q50-mu0);
	sd(idx,6) = std(q16 - (-s+mu0));
%	sd(idx,7) = std(q84 - (-s+mu0));
%	sd(idx,7) = std(q05 - (-2*s+mu0));
	sd(idx,7) = std((q84-q16)/2-s);
	sd(idx,8) = std(q16-2*q50+q84);
end

namedfigure(1,'Bias');
clf();
%for idx=1:4
%	subplot(2,2,idx);
plot(dh,[dh.^2,abs(bias)]);
set(gca,'xtick',dh);
set(gca,'xscale','log');
set(gca,'yscale','log');
legend('h^2','mean','var','skewness','kurtosis','me','16','(84-16)/2','q16-2*q50+q84');

namedfigure(2,'Standard deviation');
clf();
plot(dh,[dh.^2,sd]);
set(gca,'xtick',dh);
set(gca,'xscale','log');
set(gca,'yscale','log');
%legend('h^2','mean','var','skewness','kurtosis','me','16','84','5');
legend('h^2','mean','var','skewness','kurtosis','me','16','(84-16)/2','q16-2*q50+q84');

bias_{jdx} = bias;
end

bias_{:}
'normalised'
dh
for idx=1:length(bias_)
	bsxfun(@times,1./dh.^2,bias_{idx})
end

