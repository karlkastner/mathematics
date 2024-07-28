function w = normalwrappedpdf_width(mu,s,p)
	if (nargin()<3)
		p = 0.5;
	end
	mu = abs(mu);
	[fc, Sc] = normalwrappedpdf_mode(mu,s);
	S0 = normalwrappedpdf(0,mu,s);
	% left
	if (S0 > p*Sc)
		w = NaN;
	else
		xl = max(mu - s, 0.5*mu);
		xl = fzeros(@(x) normalwrappedpdf(x,mu,s)- 0.5*Sc,xl);
		xr = mu+s;	
		xl = fzeros(@(x) normalwrappedpdf(x,mu,s) - 0.5*Sc,xr);
	end
	w = xr - xl;
end

