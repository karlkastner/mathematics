% 2024-06-29 15:31:54.785646609 +0200
function [w,xl,xr,istruncated] = normalmirroredpdf_width(mu,s,p)
	if (nargin()<3)
		p = 0.5;
	end
	mu = abs(mu);
	[fc, Sc] = normalmirroredpdf_mode(mu,s);
	S0 = normalmirroredpdf(0,mu,s);
	% left
	if (S0 > p*Sc)
		istruncated = true;
		xl = 0;
	else
		istruncated = false;
		xl = max(mu - s, 0.5*mu);
		xl = fzero(@(x) normalmirroredpdf(x,mu,s)- 0.5*Sc,xl);
	end
	xr = mu+s;	
	xr = fzero(@(x) normalmirroredpdf(x,mu,s) - 0.5*Sc,xr);
	w = xr - xl;
end

