function x = generalized_norminv(c,mu,a,b)
	x = sign(c - 0.5).*(a.^b.*gaminv(2*abs(c-0.5),1./b)).^(1./b) + mu;
end

