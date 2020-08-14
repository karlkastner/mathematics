% Fri 13 Mar 10:42:30 +08 2020
function c = gaussfit_quantile(p,q)
	c = [mean(q),std(q)];
	c = lsqnonlin(@(c) norminv(p,c(1),c(2))-q, c);
end

