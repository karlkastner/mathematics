% 2016-02-29 13:36:34.401870431 +0100
function se = se_sample_quantiles(p,n)
	if (nargin()<2)
		n = 1;
	end
	c=normpdf(norminv(p));
	sd = 1./c.*sqrt(p.*(1-p)./n);
end

