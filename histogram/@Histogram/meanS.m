function mu = meanS(h,edge,varargin)
	mu = Histogram.momentS(h,edge,1,varargin{:});
end
