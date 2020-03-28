% 2016-03-10 16:19:07.495425422 +0100

function m = momentS(h,edge,order,method)
	if (nargin() < 4)
		method = Histogram.METHOD;
	end
	switch (method)
	case {'constant'} % constant pdf over bin
		m = constant(h,edge,order);
	case {'midpoint'}   % mid point rule, not strictly positive
		centre = Histogram.centreS(edge);
		%C = centre(1:end-1);
		C = centre;
		%w = diff(edge);
		m = sum(bsxfun(@times,h,(C.^order)),2);
	case {'gaussian'}
		m = gaussian(h,edge,order);
	case {'trapezoidal'}
		m = trapezoidal(h,edge,order);
	case {'twopoint'}
		m = twopoint(h,edge,order,2);
	case {'simpson'}
		m = simpson(h,edge,order);
%		L = edge(1:end-1);
%		C = Histogram.centreS(edge);
%		R = edge(2:end);
%		p = order;
%		m = 1/6*sum(bsxfun(@times,h,(R.^p + 4*C.^p + L.^p)),2);
%		m = 1/3*sum(bsxfun(@times,h,-0.5*L.^p +4*C.^p + -0.5*R.^p),2);	 % 4*mp - 1*trap, not positive
%		m =     2*sum(bsxfun(@times,h,iconst(L,R,order)),2) ... % 2*const - trapezoidal, not positive
%		    - 0.5*sum(bsxfun(@times,h,(R.^p + L.^p)),2);
%		m =     -sum(bsxfun(@times,h,iconst(L,R,order)),2) ... % - const + 2*mp, not positive
%		    + 2*sum(bsxfun(@times,h,(C.^p)),2);
	case {'mixture'}
		m = mixture(h,edge,order);
%	case {'linear'}
%		m = linear(h,edge,order);
%		m = 3*constant(h,edge,order)-2*linear(h,edge,order);	
	otherwise
		error('not yet implemented');
	end % switch
end % momentS

function m = simpson(h,edge,p)
	n    = size(h);
	h    = [zeros(n(1),1),h,zeros(n(1),1)];
	edge = [2*edge(1)-edge(2),edge,2*edge(end)-edge(end-1)];
	C = 0.5*(edge(2:end)+edge(1:end-1));
	width  = (edge(2:end)-edge(1:end-1));
	h = bsxfun(@times,h,1./width);
	m    = zeros(n(1),1);
	for idx=2:n(2)+1
		m = m + 1/6*(width(idx)*(C(idx-1).^p*h(:,idx-1)+2*C(idx).^p*h(:,idx)) ...
			   + width(idx+1)*(2*C(idx).^p*h(:,idx)+C(idx+1).^p*h(:,idx+1)));
	end
end

function m = linear(h,edge,p)
	n    = size(h);
	h    = [zeros(n(1),1),h,zeros(n(1),1)];
	edge = [2*edge(1)-edge(2),edge,2*edge(end)-edge(end-1)];
	centre = 0.5*(edge(2:end)+edge(1:end-1));
	width  = (edge(2:end)-edge(1:end-1));
	h = bsxfun(@times,h,1./width);
	m    = zeros(n(1),1);
	for idx=1:n(2)+1
		L = centre(idx);
		R = centre(idx+1);
		c = (h(:,idx+1)-h(:,idx))./(R-L);
		h0 = h(:,idx)-c*L; 
		m = m + 1/(p+1)*(R.^(p+1)-L^(p+1))*(h0) ... %(:,idx)-c*L) ...
		      + c/(p+2)*(R.^(p+2)-L.^(p+2));
	end
end

% histogram moment based on mixture distribution
function m = mixture(h,edges,order)
	mu    = 0.5*(edges(1:end-1)+edges(2:end));
%	sigma = 0.5*(edges(2:end)-edges(1:end-1));
	sigma = 0.5*(edges(2:end)-edges(1:end-1));
	m     = normmoment(mu,sigma,order);
	m     = sum(bsxfun(@times,h,m),2);
if (0)
	jdx = order;
	% number of histogram bins
	n = size(h);
	% width of histogram bins
	width = edges(2:end)-edges(1:end-1);
	% centre of bins
	centre = 0.5*(edges(1:end-1)+edges(2:end));
	% mean of histogram
	mu = sum(bsxfun(@times,h,mui),2);
	% central moments for normally approximated bins
	M = zeros(4,n(2));
	M(1,:) = 0;
	M(2,:) = (0.5*width).^2;
	M(3,:) = 0;
	M(4,:) = 3;
	
	m = zeros(n(1),1);
	for idx=1:n(2)
		for kdx=0:jdx
			m = m + (nchoosek(jdx,kdx)*M(kdx+1,idx))*h(:,idx).*(centre(idx) - mu).^(jdx-kdx);
		end
	end
end
end

function m = twopoint(h,edge,p,q)
	L_ = edge(1:end-1);
	R_ = edge(2:end);
	L  = q*L_ + (1-q)*R_;
	R  = (1-q)*L_ + q*R_;
	m = 0.5*sum(bsxfun(@times,h,(R.^p + L.^p)),2);
	
end

function m = trapezoidal(h,edge,p)
	m = twopoint(h,edge,p,1);
%	L = edge(1:end-1);
%	R = edge(2:end);
%	m = 0.5*sum(bsxfun(@times,h,(R.^p + L.^p)),2);
%	h_ = (R(1)-L(1));
%	m = m - 1/8*h_.^2;
end

function m = gaussian(h,edge,p)
	m = twopoint(h,edge,p,0.5-sqrt(1/2));
%	L_ = edge(1:end-1);
%	R_ = edge(2:end);
%	a = 0.5-sqrt(1/12);
%	L = a*L_ + (1-a)*R_;
%	R = (1-a)*L_ + a*R_;
%	m = 0.5*sum(bsxfun(@times,h,(R.^p + L.^p)),2);
end

function m = constant(h,edge,p)
	L = edge(1:end-1);
	R = edge(2:end);
%function int = iconst(L,R,p)
	p1 = p+1;
	int = 1./(p1*(R-L)).*( ...
		   (min(R,0).^p1 - min(L,0).^p1) ...
		 + (max(R,0).^p1 - max(L,0).^p1) ...
		);
	m = sum(bsxfun(@times,h,int),2);
%		  (min(L,0).^p1 - L.^p1) ... % from left to zero
%		+ (R.^p1 - min(L,0).^p1) ... % from zero to right
end % iconst

