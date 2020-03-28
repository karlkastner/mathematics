function s = stdS(h,edge,varargin)
	s2 = Histogram.varS(h,edge,varargin{:});
	s  = sqrt(s2);
end

