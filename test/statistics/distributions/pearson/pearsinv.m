% 2024-05-24 15:32:45.481448359 +0200
function x = pearsinv(c,mu,sd,sk,ku)
	x0 = norminv(c,mu,sd);
	x  = lsqnonlin(@(x) pearscdf(x,mu,sd,sk,ku) - c, x0);
end

