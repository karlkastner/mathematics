% 2021-09-29 15:11:53.670339080 +0200
function sd = lognpdf_std(lmu,lsd)
	%mu = lognpdf_mean(lmu,lsd);
	sd = sqrt(exp(lsd.^2)-1).*exp(lmu+0.5*lsd.^2);
end
