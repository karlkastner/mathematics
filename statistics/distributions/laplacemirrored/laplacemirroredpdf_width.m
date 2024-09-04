% 2024-06-30 11:29:29.248514283 +0200
function [w,xl,xr,istruncated] = laplacemirroredpdf_width(mu,s,p)
	if (nargin()<3)
		p = 0.5;
	end
	mu = abs(mu);
	[fc, Sc] = laplacemirroredpdf_mode(mu,s);
	S0 = laplacemirroredpdf(0,mu,s);
	% left
	if (S0 > p*Sc)
		istruncated = true;
		xl = 0;
	else
		istruncated = false;
		xl = max(mu - s, 0.5*mu);
		xl = fzero(@(x) laplacemirroredpdf(x,mu,s)- 0.5*Sc,xl);
	end
	xr = mu+s;	
	xr = fzero(@(x) laplacemirroredpdf(x,mu,s) - 0.5*Sc,xr);
	w = xr - xl;
end

