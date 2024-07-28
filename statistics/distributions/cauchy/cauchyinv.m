% 2024-05-24 16:50:35.218926551 +0200
function x = cauchyinv(c,mu,gamma)
	x = mu + gamma*tan(pi*(c-0.5));
end
