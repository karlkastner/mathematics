function mo = modeS(h,edge);
%	H      = Histogram.cdfS(h);
	centre = Histogram.centreS(edge);
	width  = abs(edge(2:end)-edge(1:end-1));

	% scale for variable bin width
	h = bsxfun(@times,h,1./width);

	% allocate memory
	n  = size(h,1);
	mo = zeros(n,1,class(h));

	% maximum
	[mv mdx] = max(h,[],2);

	% to avoid out of bound errors empty cells are padded to the ends
	h = [zeros(n,1), h, zeros(n,1)];

	% extrapolate exterior cell centre
	centre = [2*centre(2)-centre(1) centre(:)' 2*centre(end)-centre(end-1)];

	% maximum shifts one to the right because of the padding
	mdx = mdx+1;

	% find maximum by approximating peak as a quadratic
	for idx=1:n
		hi = h(idx,mdx(idx)-1:mdx(idx)+1);
		hi = hi(:);
		x = centre(mdx(idx)-1:mdx(idx)+1);
		x = x(:);
		A = [ones(3,1) x x.*x];
		c = A \ hi;
		x0 = -0.5*c(2)/c(3);
		mo(idx) = x0;
	end
end

