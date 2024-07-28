function [w,xl,xr,istruncated] = cauchywrappedpdf_width(mu,s,p)
	if (nargin()<3)
		p = 0.5;
	end
	mu = abs(mu);
	[fc, Sc] = cauchywrappedpdf_mode(mu,s);
	S0 = cauchywrappedpdf(0,mu,s);
	% left
	if (S0 > p*Sc)
		istruncated = true;
		xl = 0;
	else
		istruncated = false;
		xl = max(mu - s, 0.5*mu);
		xl = fzero(@(x) cauchywrappedpdf(x,mu,s)- 0.5*Sc,xl);
	end
	xr = mu+s;	
	xr = fzero(@(x) cauchywrappedpdf(x,mu,s) - 0.5*Sc,xr);
	w = xr - xl;
end

