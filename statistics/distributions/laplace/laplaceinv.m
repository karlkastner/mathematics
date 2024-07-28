% 2024-07-08 13:08:30.625946361 +0200
function x = laplaceinv(q,mu,s)
	flag = (q<1/2);
	x = mu + s*(log(2*q).*flag - log(2-2*q).*(~flag));
end

