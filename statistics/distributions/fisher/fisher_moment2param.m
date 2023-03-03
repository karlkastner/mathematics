% Mon 19 Dec 21:16:09 CET 2022
function [d1,d2] = fisher_moment2param(mu,sd)
	d2 = 2*mu/(mu-1);
	d1 = (d2 - 2)/((sd^2*(d2 - 2)^2*(d2 - 4))/(2*d2^2) - 1);
end

