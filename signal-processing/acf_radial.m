% 2021-06-23 16:22:53.456353906 +0200
function mu = r_acf(f)
	[mu,se,r,n] = r_spectrum(f);
	mu = mu/mu(1);
end
