% Wed 15 May 21:28:49 CEST 2024
% Karl Kastner, Berlin
% TODO this assumes equally wide histogram bins
function par = sn2_fit(x,h,par0)
	x = double(x);
	h = double(h);
	x = cvec(x);
	h = cvec(h);
	h = h/sum(h);
	if (nargin()<3)
	% initial guess
	n1 = find(cumsum(h)>0.5,1,'first')-1;
	fdx = (1:n1);
	q = quantile(h,[0.1,0.9]);
	p1  = sum(h(fdx));
	mu1 = sum(h(fdx).*x(fdx))/p1;
	sd1 = sqrt(sum(h(fdx).*(x(fdx)-mu1).^2)/p1);
	sk1 = 0;
	fdx = (n1+1:length(h));
	p2 = 1-p1;
	mu2 = sum(h(fdx).*x(fdx))/p2;
	sd2 = sqrt(sum(h(fdx).*(x(fdx)-mu2).^2)/p2);
	sk2 = 0;
	par0 = [p1,mu1,mu2,sd1,sd2,sk1,sk2];
	else
	par0 = double(par0);
	end
	smax = 0.995;
	lb  = [0,-inf,-inf,0,0,-smax,-smax];
	ub  = [inf,inf,inf,inf,inf,smax,smax];
	h = h/sum(h*(x(2)-x(1)));
	par = lsqnonlin(@(p) skew_normal_mixture_pdf(x,[p(1),1-p(1)],p(2:3),p(4:5),p(6:7)) - h,par0,lb,ub);
end % sn2_fit

