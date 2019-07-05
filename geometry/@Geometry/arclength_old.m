% Thu Apr 14 16:31:37 CEST 2016
% Karl Kastner, Berlin
%
%% arc length of a two dimensional function
% TODO deprecated, use arc length
function s = integrate_path(x,y,order)
	x = rvec(x);
	y = rvec(y);

	switch (order)
%	case {1} % chord length
%		s = [0,cumsum(hypot(diff(x),diff(y)))];
	case {12}
		% padd to uneven number of elements
		flag = mod(length(x),2);
		if (0 == flag)
			x(end+1) = x(end);
			y(end+1) = y(end);
		end
%		s   = [0,cumsum(hypot(diff(x),diff(y)))];
%		s1  = s(1:2:end);
%		s2  = [0;cumsum(hypot(diff(x(1:2:end)),diff(y(1:2:end))))];
%		s12 = 1/3*(4*s1-s2);
		% there are now only 
%		s = reshape(s,
		
		% segment length
		dS1 = hypot(diff(x),diff(y));
		dS1 = [dS1(1:2:end);dS1(2:2:end)];
		% segment length with skipping one point
		dS2 = hypot(diff(x(1:2:end)),diff(y(1:2:end)));
		% unequal step-width factor
		a = dS1(1,:)./dS1(2,:);
		a = 1;
		%f = (1+a.^2)./(1+a).^2;
		% higher order accurate segment length
		dS1_ = sum(dS1);
		dS12 = (2*dS1_./(1+a.^2)-dS2./(1+a).^2)./(2./(1+a.^2)-1./(1+a).^2);;


		%dS12 = 1/3*(4*dS1_-dS2);
		%dS12 = (2*dS1_-dS2);
		% scale up
		dS  = bsxfun(@times,dS1,dS12./dS1_);
		% sum up pieces
		s = [0,cumsum(dS(:)')];
		if (0 == flag)
			s = s(1:end-1);
		end
%		dS = reshape(dS,
%		X = [x(1:2:end); x(2:2:end); x(3:3:end)];
%		Y = [y(1:2:end); y(2:2:end); y(3:3:end)];
		% close tail
	case {2} % quadratic
		l = length(x);
		m = 2*floor((l-1)/2);
		dS = zeros(size(x));
		dS(2:m+1) = int3(x(1:m+1),y(1:m+1));
		dS(end-1:end) = int3(x(end-2:end),y(end-2:end));
		s = cumsum(dS);
		if (length(x) ~= length(s))
			error('here')
		end
	otherwise % chord length (linear (1))
		s = [0, cumsum(hypot(diff(x),diff(y)))];
	end
end

% as a circle
	% find circumcentre and radisu of three points
	% get angle for A,B, centre
	% get arc length as R*alpha_rad
	% get angle for B,C, centre
	% get arc length as R*beta_rad

% TODO, if central point not convex, fall back to lower order scheme
%if (dx > dy)
%	make y a function of x and integrate
%else
%	make x a function of y and integrage

	

function dS = int3(x,y)
	xl = x(1:2:end-2);
	xc = x(2:2:end-1);
	xr = x(3:2:end);
	yl = y(1:2:end-2);
	yc = y(2:2:end-1);
	yr = y(3:2:end);

% x = c(1) + c(2)*t + c(3)*t^2
cx = [xc; xr/2 - xl/2; xl/2 - xc + xr/2];
cy = [yc; yr/2 - yl/2; yl/2 - yc + yr/2];

% dx^2 + dy^2 = a + b*t + c*t^2
a = cx(2,:).^2 + cy(2,:).^2;
b = 2*cx(2,:).*cx(3,:) + 2*cy(2,:).*cy(3,:);
c = cx(3,:).*cx(3,:) + cy(3,:).*cy(3,:);

I = [(a.^(1./2).*b)./(4.*c) - (b./(4.*c) - 1./2).*(a - b + c).^(1./2) ...
   + (log(b./(2.*c.^(1./2)) + a.^(1./2)).*(- b.^2./4 + a.*c))./(2.*c.^(3./2)) ... 
   - (log((b./2 - c)./c.^(1./2) + (a - b + c).^(1./2)).*(- b.^2./4 + a.*c))./(2.*c.^(3./2));
     (b./(4.*c) + 1./2).*(a + b + c).^(1./2) ...
   + (log((b./2 + c)./c.^(1./2) + (a + b + c).^(1./2)).*(- b.^2./4 + a.*c))./(2.*c.^(3./2)) ...
   - (a.^(1./2).*b)./(4.*c) ...
  - (log(b./(2.*c.^(1./2)) + a.^(1./2)).*(- b.^2./4 + a.*c))./(2.*c.^(3./2))];
 	dS = reshape(I,1,[]);
	size(dS)
	size(xl)
end

