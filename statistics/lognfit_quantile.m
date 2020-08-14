% Fri 13 Mar 10:42:30 +08 2020
% function c = lognfit_quantile(p,q)
function c = lognfit_quantile(p,q)
	c = [mean(log(q)),std(log(q))];
	opt=  optimset('display','off');
	c = lsqnonlin(@(c) logninv(p,c(1),c(2))-q,c,[],[],opt);
end
