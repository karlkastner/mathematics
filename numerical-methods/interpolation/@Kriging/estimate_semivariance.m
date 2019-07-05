% Sun Jun 15 23:57:29 WIB 2014
%% estimate the parameter of the semivariance model for Kriging interpolation
function obj = estimate_semivariance(obj, Xs, Vs, XY_)
	% TODO, make this a class member
	ncalibrate = 1000;
	obj.Rmax = 1500;

	% filter invalid and duplicate values (these may crash the tree)
	% TODO filtration should be done after construction
	fdx = isfinite(prod(Xs,2));
	Xs  = Xs(fdx,:);
	Vs  = Vs(fdx,:);
	[Xs udx] = unique(Xs,'rows');
	Vs = Vs(udx);	

	% efficient estimate of the semivariance model
	% TODO, this is only implemented for the exponential model so far
	% TODO, this is only implemented for 2d so far
	% estimate the first parameter
	% by drawing random pairs of samples
	% here as many samples as points are drawn
	% note : duplicate samples do not affect the parameters, but their covariance
	n = length(Xs);
	idx = randi(n,n,1);
	jdx = randi(n,n,1);
	% sampled semivariance and its variance :)
	sv = (Vs(idx) - Vs(jdx)).^2;
	% remove outliers
	sdx = find(sv<quantile(sv,0.95));
	s2 = 0.5*mean(sv(sdx));
	% uncertainty of the variance parameter
	sparam(1) = 0.25/length(sdx)*var(sv(sdx));
	% estimate the exponential parameters from the nearest neighbours
	% for each point of the calibration subset
	R2max = obj.Rmax*obj.Rmax;
	b = [];
	A = [];
	c =[];
	% choose a subset of points for calibration
	jd = randi(n,ncalibrate,1);
	for idx=1:ncalibrate
		% TODO use qtree to speed this up
		d = obj.dist(Xs,Xs(jd(idx),:));
		d2 = d(:,1).^2 + d(:,2).^2;
		% exclude itself
		id = find(d2 < R2max & (1:n)' ~= jd(idx));
		b_ = 0.5*log((1-0.5/s2*(Vs(id)-Vs(jd(idx))).^2).^2);
		% switch back to euclidean distance
		A_ = [-abs(Xs(id,1)-Xs(jd(idx),1)) -abs(Xs(id,2)-Xs(jd(idx),2))];
		b = [b;b_];
		A = [A;A_];
		c_ = 0.5/s2*(Vs(id)-Vs(jd(idx))).^2;
		c = [c;c_];
	end
	% regress the parameter
	p = A\b;
	% exclude outliers
	res = A*p - b;
	res2 = res.*res;
	fdx = find(res2 < quantile(res2,0.95));
	% regress parameters again
	p = A(fdx,:)\b(fdx);
	err = A*p-b;
	param = [s2 p'];

%		qtree_ = Qtree_(Xs(:,1),Xs(:,2));
%		[ndx dmin] = qtree_.nearest_neighbour(Xs(:,1),Xs(:,2),2);
%		ndx = ndx(:,2);
%		dmin = dmin(:,2);
%		% only use the lower 95%
%		q = quantile(dmin,0.95);
%		fdx = find(dmin < q);
%%		% set up the regression matrix and solve for parameters
%		A = [abs(Xs(fdx,1) - Xs(ndx(fdx),1)) abs(Xs(fdx,2) - Xs(ndx(fdx),2))];
%		b = (Vs(fdx) - Vs(ndx(fdx))).^2/(2*s2);
%		param(2:3,1) = A \ b;
%		err = A*param(2:end) - b;
	covp = err'*err/(length(fdx)-2)*inv(A(fdx,:)'*A(fdx,:));
	sparam(2:3,1) = diag(covp);
	obj.param  = param;
	obj.sparam = sqrt(sparam);

	binx = linspace(0,obj.Rmax*obj.s(1),12.5);
	biny = linspace(0,obj.Rmax*obj.s(2),12.5);
	px  = bin1d(-A(:,1),c,binx);
	py = bin1d(-A(:,2),c,biny);
	h=bin2d(-A(:,1),-A(:,2),c,binx,biny);
	%linspace(0,500,10),linspace(0,500,10));
	clf();
	subplot(2,1,1)
	plot(binx(1:end-1),px);
	hold on
	plot(biny(1:end-1),py,'g');
	subplot(2,1,2)
	bar3(h);
%plot([h(:,1) h(1,:)'])
%	par = A \ log(1-2*b)
%'molch'
%	pause
end % estimate_semivariance

