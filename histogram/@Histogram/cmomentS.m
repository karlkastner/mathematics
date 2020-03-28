function m = cmomentS(h,edge,order,varargin);
	m = 0;
	n  = order;
	m1 = Histogram.momentS(h,edge,1,varargin{:});
	for idx=0:n
		m = m + (-1)^(n-idx)*nchoosek(n,idx)*Histogram.momentS(h,edge,idx,varargin{:}).*m1.^(n-idx);
	end
end % cmomentS

