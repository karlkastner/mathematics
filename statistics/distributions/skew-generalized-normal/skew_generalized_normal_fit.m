% Wed 25 Nov 13:35:34 +08 2020
% Karl Kastner, Berlin
function [c,pdf] = skew_generalized_normal_fit(h,x,c,w)
	if (nargin()<3)
		mu = rvec(h)*cvec(x)
		sd = sqrt(rvec(h)*(cvec(x)-mu).^2);
		c = [mu,sd,-1];
	end
	if (nargin()<4)
		w = 1;
	end
	mynorm = @(x) x/sum(x);
	%pdf    = @(x,c) mynorm(skew_generalized_normpdf(x,c(1),c(2),tan(c(3)/(pi*2)),c(4))); %sign(c(3))*10.^abs(c(3)),c(4)));
	%pdf    = @(x,c) mynorm(skew_generalized_normpdf(x,c(1),c(2),c(3)/1e4,c(4))); %sign(c(3))*10.^abs(c(3)),c(4)));
	pdf    = @(x,c) mynorm(skew_generalized_normpdf(x,c(1),c(2),c(3),c(4))); %c(4))); %sign(c(3))*10.^abs(c(3)),c(4)));
	opt = struct();
	opt.Algorithm = 'levenberg-marquardt';
	lbg = [];
	ubg = [];
	c     = lsqnonlin(@(c) w.*(h - pdf(x,[c,0])),c,lbg,ubg,opt);
	c     = lsqnonlin(@(c) w.*(h - pdf(x,c)),[c,1],lbg,ubg,opt);
end
