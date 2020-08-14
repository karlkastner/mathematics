% 2020-02-26 10:40:14.864663112 +0800
function param = gaussfit_quantile(p,q,c0)
	if (nargin()<3)
	c0    = [mean(q),std(q)]; 
	end
	c = c0;
	opt = optimset('Display','off');
	[param,rn,res,flag] = lsqnonlin(@(c) norminv(p,c(1),c(2)) - q,c0,[],[],opt);
end
