% 2015-02-06 14:58:01.684729867 +0100
% Karl Kastner, Berlin
%% fit slope and intercept to a set of sample with the Theil-Sen method
%% 
%% c     : confidence interval c = 2*ns*normcdf(1) for ns-sigma intervals
%% param : itercept and slope
%% P : confidence interval
function [param, res, obj] = fit(obj,X,Y,W,P)
	if (nargin()<4)
		W = [];
	end
	% TODO, put inference into predict
	if (nargin()<5)
		P = 0.95;
	end

	if (isvector(X))
		X = cvec(X);
	end
	if (isvector(Y))
		Y = cvec(Y);
	end

	% number of samples
	n = length(X);
	if (n < 2)
		param = [NaN;NaN];
		return;
	end

	% number of sample pairs
	np    = 0.5*n*(n-1);
	% critical value
	% two sided band, so devide by 2
	%z     = norminv(0.5*c);
	z = norminv(1-0.5*(1-P));

	%
	% constant term
	%

	% ranks for the median confidence intervals
	% TODO the odd even thing is not in the original formula,
	% but makes sense in my opinion
%	d = z*sqrt(n);
%	if (1 == mod(n,2))
%		rankl = 0.5*(n - d);
%		ranku = 0.5*(n + d) + 1;
%	else
%		rankl = 0.5*(n - d - 1);
%		ranku = 0.5*(n + d + 1);
%	end
%	% limitation is necessary if there are insufficient many samples
%	pl = max(0,rankl/n);
%	pm = 0.5;
%	pu = min(1,ranku/n);
%
%	qx    = quantile(X,[pl pm pu]);
%	qy    = quantile(Y,[pl pm pu]);
	[m, void, l, u] = median_man([X, Y],P);
	qx = [l(1), m(1), u(1)];
	%[m void l u] = median_man([X Y],P);
	qy = [l(2:end); m(2:end); u(2:end)];

	%
	% slope
	% 

	% ranks for slope confidence intervals
	% see helsel and hirsch (2002)
	dp     = z*sqrt(n*(n-1)*(2*n+5)/18);
	if (mod(np,1))
		rankl = 0.5*(np - dp);
		ranku = 0.5*(np + dp) + 1;
	else
		rankl = 0.5*(np - dp - 1);
		ranku = 0.5*(np + dp + 1);
	end
	% limitation is necessary if there are insufficient many samples
	pl = max(0,rankl/np);
	pm = 0.5;
	pu = min(1,ranku/np);

	[slope, qs] = obj.slope(X,Y,[],[pl, pm, pu]);
%	XX   = repmat(X,1,n);
%	YY   = repmat(Y,1,n);
%	dx   = XX-XX';
%	dy   = YY-YY';
%	flag = triu(true(n),1);
%	dy_dx = dy(flag(:))./dx(flag(:));
%	qs    = quantile(dy_dx,[pl pm pu]);


	% TODO, why not median of the interectps of all pairs?

	% Note that there exists also an alternative definition,
	% which is median(y - slope*x) that does not exactly yield the same result
	% note that the current definition of the intercept is recommended by
	% some authors, but has obviously convergence problems if sampling range is increased
	% TODO use hodges lehman (!)
%	slope = qs(2);
	mey = qy(2,:);
	mex = qx(2);
	% intercept
	intercept = mey - slope*mex;
%	qi = quantile(Y - slope*X,
	% order of X*slope important for matrix Y
	[m, void, l, u] = median_man(Y - X*slope, P);
	qi = [l, m, u];

	param      = [intercept; slope];
	obj.param  = param;
	obj.params = 0.5*[qi(3)-qi(1),qs(3)-qs(1)];
	obj.qx    = qx;
	obj.qy    = qy;
	obj.qs    = qs;
	obj.qi    = qi;

	% TODO this is OLS not TS	
	Yp  = obj.predict(X);
	res = Yp - Y;
	n = length(Y);
	obj.serr = sqrt(n/(n-2))*rms(res);


	r2 = struct();
	% coefficient of correlation
	r2.pearson  = 1 - obj.serr.^2./var(Y);
	% identical to:
	% r2_ = corr(Y,Yp,'type','Pearson')^2;
	% r2.pearson = 1 - (ne-1)/(ne-np)*(1-r2_)
	if (obj.extended_statistics)
		np = 2;
		ne = length(Y);
		r2_ = spearman_to_pearson(corr(Y,Yp,'type','Spearman'))^2;
		r2.spearman = (ne-1)/(ne-np)*(1-r2_);
		r2_ = kendall_to_pearson(corr(Y,Yp,'type','Kendall'))^2;
		r2.kendall  = (ne-1)/(ne-np)*(1-r2_);
		r2_ = hodges_lehmann_correlation(Y,Yp);
		r2.hodges_lehmann  = (ne-1)/(ne-np)*(1-r2_);
		obj.mad     = median(abs(res));
	end

	obj.r2 = r2;
	%R2 = 1 - obj.serr.^2./(var(Y));

end % regress

