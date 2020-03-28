% So 13. Mär 12:46:47 CET 2016
function s2 = varS(h,edge,varargin);
%	s2 = Histogram.
	s2 = Histogram.cmomentS(h,edge,2,varargin{:});
%	mu = Histogram.meanS(h,edge,varargin{:});
%	x2 = Histogram.momentS(h,edge,2,varargin{:});
%	s2 = x2 - mu.^2;
%	mu  = histmean(h,centre);
%	d   = repmat(centre(:)',size(h,1),1) - repmat(mu,1,size(h,2));
%	s2  = sum(h.*(d.*d),2);

	if(0) %nargin() > 2)
	switch(varargin{1})
	case {'midpoint'}
%	if (nargin > 2 && strcmp(varargin,'trapezoidal'))
		s2 =s2 - 1/12*(edge(2)-edge(1))^2;
	case {'trapezoidal'}
%	if (nargin > 2 && strcmp(varargin,'trapezoidal'))
		s2 =s2 - 1/3*(edge(2)-edge(1))^2;
	case {'constant'}
		s2 = s2 - 1/6*(edge(2)-edge(1)).^2;
	case {'gaussian'}
		s2 = s2 - 1/6*(edge(2)-edge(1)).^2;
	case {'simpson'}
		s2 = s2 - 1/12*(edge(2)-edge(1)).^2;
	end


	% bias correction
	% TODO, this assumes constant bin width
	% TODO, this should go into the moment computation
	if (Histogram.SHEPPARD)
		centre = Histogram.centre;
		dh  = centre(2)-centre(1);
		s2  = max(sqrt(eps),s2 - 1/12*dh.^2);
	end
end

