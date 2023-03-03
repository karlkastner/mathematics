% 2021-09-29 15:11:48.878175235 +0200
function mu = logn_mean(lmu,lsd)
	mu = exp(lmu + 0.5*lsd.^2);
end


