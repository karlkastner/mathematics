function sk = skewnessS(h,edge,varargin)
	m3 = Histogram.cmomentS(h,edge,3,varargin{:});
	s  = Histogram.stdS(h,edge,varargin{:});
	% even moments are unbiased, no Sheppard correction necessary
	sk = m3./(s.^3);	

%	mu  = histmean(h,centre);
%	s   = histstd(h,centre);
%	d   = repmat(centre(:)',size(h,1),1) - repmat(mu,1,size(h,2));
%	mu3 = sum(h.*(d.*d.*d),2);
%	sk  = mu3./(s.*s.*s);
end

