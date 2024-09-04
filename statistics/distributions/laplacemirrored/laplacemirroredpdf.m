% Sat 29 Jun 10:29:30 CEST 2024
function p = laplacemirroredpdf(x,mu,s)
	p = 0.5*(laplacepdf(x,mu,s) + laplacepdf(x,-mu,s));
end

