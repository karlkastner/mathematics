function S = gammirroredpdf(x,a,b)
	S = 0.5*(gampdf(x,a,b) + gampdf(-x,a,b));
end
