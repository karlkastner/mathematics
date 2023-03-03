% 2021-09-29 15:11:53.670339080 +0200
function sd = logn_std(lmu,lsd)
	mu = logn_mean(lmu,lsd);
	sd = sqrt(exp(lsd.^2)-1).*exp(lmu+0.5*lsd.^2);
end
