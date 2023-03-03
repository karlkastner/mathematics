
function h = lognpdf_entropy(lmu,lsd)
	h = log2(sqrt(2*pi)*lsd*exp(lmu+0.5));
end
