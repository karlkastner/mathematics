% Mon 13 Feb 13:42:27 CET 2023
function [fc,Sc] = laplacemirroredpdf_mode(f0,s)
	%dS_df = - exp(-abs(f + f0)/s) - (exp(-abs(f - f0)/s)*sign(f - f0))
	dS_df = @(f) - exp(-abs(f - f0)/abs(s))*sign(f - f0) - (exp(-abs(f + f0)/abs(s)).*sign(f + f0));
	fc = fzero(dS_df,f0);
	Sc = laplacemirroredpdf(fc,f0,s);
end

