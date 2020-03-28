% Fri 13 Mar 10:42:30 +08 2020
function c = lognfitpq(p,q)
	c = [mean(q),1];
	c = lsqnonlin(@(c) logninv(p,c(1),c(2))-q,c);
	%c = lsqnonlin(@(c) logncdf(q,c(1),c(2))-p,c);
end
