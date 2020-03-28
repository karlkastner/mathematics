function kurt = kurtosisS(h,edge,varargin)
	m4 = Histogram.cmomentS(h,edge,4,varargin{:});
	s2 = Histogram.varS(h,edge,varargin{:});

	% apply Sheppard's bias correction
	% TODO, this should go into moment function
	if (Histogram.SHEPPARD)
		centre = Histogram.centreS(edge);
		dh  = centre(2)-centre(1);
		% should this be s2 or m2 ?
		m2 = s2;
		%m4  = m4 - 0.5*m2*dh^2 + 7/240*dh^4;
	end
	kurt = m4./(s2.*s2);
	% empirical
	%kurt = kurt - 0.5*dh^2./s2 - 1/80*dh^4./s2.^2;

%	mu1 = histmean(h,centre);
%	s2  = histvar(h,centre);
%	d1  = repmat(centre(:)',size(h,1),1) - repmat(mu,1,size(h,2));
%	d2  = d1.*d1;
%	d4  = d2.*d2; 
	% TODO, this assumes constant bin-width
%	mu2 = sum(h.*d2,2);
%	mu4 = sum(h.*d4,2);
end

