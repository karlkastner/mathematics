function [a,b] = lognmirroredpdf_mode2par(xc,Sc)
	[a,b] = lognpdf_mode2par(xc,2*Sc); 
end

