% Sat 10 Dec 15:44:16 CET 2016
% Karl Kastner, Berlin
%% smooth a parametric function given in x-y coordinates
function [x y] = smooth_1d(x,y,order)
	if (nargin() < 3)
		order = 1;
	end
	method = 0;
	method = 'quintic';
	method = 'ip';
	im = 'spline';
	switch (method)
	case {'ip'}
		% TODO this requires better integration with the simpson rule
		% expand x and y locally as functions in T
		%dx = diff(x);
		%dy = diff(y);
		%ds = hypot(dx,dy);
		%S  = [0;cumsum(ds)];
		S = Geometry.arclength(x,y);
		Si = S;
		for idx=1:order
			Si(2:end-1) = 0.5*(Si(1:end-2)+Si(3:end));
		end
		x = interp1(S,x,Si,im);
		y = interp1(S,y,Si,im);
	case {-1}
		% at mid point with orthogonal displacement such that path length remains the same
		% ill conditioned
	case (0) % regress a circle and put to centre
		% preserves radius, but leads to inversion in sharp bends
		% and if ds varies in bends, sharp edges
		x = cvec(x);
		y = cvec(y);
		[x0 y0 R1] = Geometry.tri_excircle([x(1:end-2),x(2:end-1),x(3:end)], ...
		                                  [y(1:end-2),y(2:end-1),y(3:end)]);
		xc = 0.5*(x(1:end-2)+x(3:end));
		yc = 0.5*(y(1:end-2)+y(3:end));
		R2 = hypot(xc-x0,yc-y0);
		xc = x0 + (xc-x0).*R1./R2;
		yc = y0 + (yc-y0).*R1./R2;
		x(2:end-1)=xc;
		y(2:end-1)=yc;
	case ('quintic')
		x = cvec(x);
		y = cvec(y);
		for idx=3:length(x)-2
			dl = hypot(diff(x(idx+(-2:2))),diff(y(idx+(-2:2))));
			l = [0;cumsum(dl)]
			l = l-l(end)/2
			A = vander_1d(l,3);
	%		Ai = inv(A);
			%cx = Ai*x(idx+(-2:2));
			%cy = Ai*y(idx+(-2:2));
			cx = A \ x(idx+(-2:2));
			cy = A \ y(idx+(-2:2));
			x_(idx)=cx(1);
			x_(idx)=cx(1);
			y_(idx)=cy(1);
		end % for idx
		x(3:end-2) = x_(3:length(x)-2);
		y(3:end-2) = y_(3:length(x)-2);
	case ('quadxy')
		% preserves length
		% makes ds equal
		% but increases noise in transversal direction
		n = length(x);
		x= rvec(x);
		y= rvec(y);
		% make x and y quadratic in chord length
		dx = diff(x);
		dy = diff(y);
		dl = hypot(dx,dy);
		A = zeros(3,3,n-2);
		A(1,1:3,:) = 1;
		A(2,1,:) = -dl(1:end-1);
		A(2,3,:) =  dl(2:end);
		A(3,1,:) =  dl(1:end-1).^2;
		A(3,3,:) =  dl(2:end).^2;
		A = transpose3x3(A);
		Ai = inv3x3(A);
		bx = [x(1:end-2);
	              x(2:end-1);
	              x(3:end)];
		cx = matvec3x3(Ai,bx);
		cx = squeeze(cx);
		by = [y(1:end-2);
	              y(2:end-1);
	              y(3:end)];
		cy = matvec3x3(Ai,by);
		cy = squeeze(cy);
		% expand
		l0 = 0.5*(dl(2:end)-dl(1:end-1));
	%	l0 = 0.*l0;
		xc = cx(1,:) + cx(2,:).*l0 + cx(3,:).*l0.^2;
		yc = cy(1,:) + cy(2,:).*l0 + cy(3,:).*l0.^2;
		x(2:end-1) = xc;
		y(2:end-1) = yc;
	case {3}	% midpoint with displacement, local
	% this smoothes ds, but not dy,dx

	% odd even
	for idx=0:1

	x = rvec(x);
	y = rvec(y);
	% put point into the centre of its neighbours, but apply a path length correction
	xl = x(idx+(1:2:end-2-idx));
	xc = x(idx+(2:2:end-1-idx));
	xr = x(idx+(3:2:end-idx));
	yl = y(idx+(1:2:end-2-idx));
	yc = y(idx+(2:2:end-1-idx));
	yr = y(idx+(3:2:end-idx));
	%xl = x(1:end-2);
	%xc = x(2:end-1);
	%xr = x(3:end);
	%yl = y(1:end-2);
	%yc = y(2:end-1);
	%yr = y(3:end);

	% length
	dxl   = xc-xl;
	dyl   = yc-yl;
	dll  = hypot(dxl,dyl);
	dxr   = xc-xr;
	dyr   = yc-yr;
	dlr  = hypot(dxr,dyr);

%	dl  = hypot(dx,dy);
	dxc = xl-xr;
	dyc = yl-yr;
	dlc = hypot(dxc,dyc);
	% base points
	xyb = Geometry.base_point([xc; yc],[xl;yl],[xr;yr]);
%	hold on
%	plot(xyb(1,:),xyb(2,:),'b.');
%	plot([xl xr],[yl yr],'b-');
	% direction from base point to in between point
	% this is orthogonal to left-right, but preserves sign
	ox = xc-xyb(1,:);
	oy = yc-xyb(2,:);
	% normalise
	h  = hypot(ox,oy);
	ox = ox./h;
	oy = oy./h;
	% centre
	xcn = 0.5*(xl+xr);
	ycn = 0.5*(yl+yr);
%	plot(xcn,ycn,'g.');
	r2 = (0.5*(dll+dlr)).^2 - (0.5*dlc).^2;
%	r = -r;
	% displacement in orthogonal direction that preserves the path length
%ox^2*r^2 - ox^2*ycn^2 + 2*ox*oy*xcn*ycn + oy^2*r^2 - oy^2*xcn^2
%r^2 - ox^2*ycn^2 + 2*ox*oy*xcn*ycn - oy^2*xcn^2
%	a = -(ox.*xcn + oy.*ycn - sqrt(r.^2 - ox.^2.*ycn.^2 + 2.*ox.*oy.*xcn.*ycn - oy.^2.*xcn.^2))./(ox.^2 + oy.^2);
%	a = -a;
	a = sqrt(r2); %./(ox^2 + oy^2)^(1/2)
	xcn = xcn + a.*ox;
	ycn = ycn + a.*oy;
	% (if this is alternated odd/even, then the path length is indeed preserved)
	%x(2:end-1) = xcn;
	%y(2:end-1) = ycn;
	x(idx+(2:2:end-1-idx)) = xcn;
	y(idx+(2:2:end-1-idx)) = ycn;

	end
%	sum(dl)
%	sum(hypot(diff(x),diff(y)))
	end
%	n = length(x);
%	sum(hypot(diff(x),diff(y)))
%
%
%
%	dx  = diff(x);
%	dy  = diff(y);
%	l   = hypot(dx,dy);
%	lc  = l(1:end-1) + l(2:end);
%	dxc = x(3:end)-x(1:end-2);
%	dyc = y(3:end)-y(1:end-2);
%	l_ = hypot(dxc,dyc);
%	dxc = dxc./l_;
%	dyc = dyc./l_;
%	% orthogonal direction
%	odx = -dyc;
%	ody =  dxc;
%	% orthogonal displacement
%	o  = 0.5*sqrt(lc.^2 - l_.^2);
%	% the sign depends on the original direction
%	R = zeros(2,2,n-2);
%	R(1,1,:) =  dxc;
%	R(1,2,:) =  dyc;
%	R(2,1,:) = -dyc;
%	R(2,2,:) =  dxc;
%	r = matvec2x2(R,[xc-xl;yc-yl]);
%	sig = -sign(r(1,:))';
%%	matvec2x2(R,[dxc;dyc])
%	
%
%	% displace to preserve path length
%	xc = xcn + sig.*odx.*o;
%	yc = ycn + sig.*ody.*o;
%
%	x(2:end-1)=xc;
%	y(2:end-1)=yc;
%
%	sum(hypot(diff(x),diff(y)))
%end

if (0)
	n = length(x);
	x = rvec(x);
	y = rvec(y);

if (0)
	% get transformation matrix
	% T = [xl xr; yl yr] / [-1 1; 0 0];
	A=zeros(2,2,n-2);
	A(1,1,:) = x(1:end-2); % xl
	A(1,2,:) = x(3:end);   % xr
	A(2,1,:) = y(1:end-2); % yl
	A(2,2,:) = y(3:end);   % yr
	B = zeros(2,2,n-2);
	B(1,1,:) = -1;
	B(1,2,:) =  1;
	T = mtimes2x2(B,inv2x2(A));

	% transform mitpoint
	t = matvec2x2(T,[x(2:end-1);y(2:end-1)]);
else
	xl = x(1:end-2);
	xc = x(2:end-1);
	xr = x(3:end);
	yl = y(1:end-2);
	yc = y(2:end-1);
	yr = y(3:end);


	% translate
	x0 = 0.5*(xl+xr);
	y0 = 0.5*(yl+yr);
	xl = xl-x0;
	xr = xr-x0;
	xc = xc-x0;
	yl = yl-y0;
	yr = yr-y0;
	yc = yc-y0;
	% scale
	s  = 1./sqrt(xl.^2 + yl.^2);
	xl = s.*xl;
	xc = s.*xc;
	xr = s.*xr;
	yl = s.*yl;
	yc = s.*yc;
	yr = s.*yr;
	% rotate
	R = zeros(2,2,n-2);
	R(1,1,:) =   xl;
	R(2,2,:) =   xl;
	R(1,2,:) =   yl;
	R(2,1,:) =  -yl;
	t = matvec2x2(R,[xc;yc]);
%	matvec2x2(R,[xl;yl])
%	matvec2x2(R,[xr;yr])
end

	% now the line is 1-dimensional, determine coefficients
	A = zeros(3,3,n-2);
%	A = [1 -1 1;
%            1  t(1) t(1)^2;
%           1  1 1];

	A(1,1,:) = 1; A(1,2,:) =     -1; A(1,3,:) = 1;
	A(2,1,:) = 1; A(2,2,:) = t(1,:); A(2,3,:) = t(1,:).^2;
	A(3,1,:) = 1; A(3,2,:) =      1; A(3,3,:) = 1;

%	A(1,1,:) =  1; A(1,2,:) =      1;    A(1,3,:) = 1;
%	A(2,1,:) = -1; A(2,2,:) = t(1,:);    A(2,3,:) = 1;
%	A(3,1,:) =  1; A(3,2,:) = t(1,:).^2; A(3,3,:) = 1;

	b = [zeros(1,n-2); t(2,:); zeros(1,n-2)];
	c = matvec3x3(inv3x3(A),b);
	
	% expand at 0 (the midpoint)
	% ty = [1 tx tx^2]*c, tx = 0
	t0 = [zeros(1,n-2); c(1,:)];
%	t0 = t;
%	t0 = [c(1,:);zeros(1,n-2)];
	
%	t  = linspace(-1,1,100);
	%t0 = c(1) + c(2)*t + c(3)*t.^3;
	% retransform
if (0)
	xy = matvec2x2(inv2x2(T),t0);
else
	% rotate
	xy = matvec2x2(inv2x2(R),t0);
	% scale
	si = 1./s;
	xy(1,:) = si.*xy(1,:);
	xy(2,:) = si.*xy(2,:);
	% translate
	xy(1,:) = xy(1,:) + x0;
	xy(2,:) = xy(2,:) + y0;
end
	
	x(2:end-1)=xy(1,:)';
	y(2:end-1)=xy(2,:)';
         
	% move bifurcation points to average of neighbours
	% TODO
end

end % smooth_1d

% Backwater for hole and riffle -> what is the depth over a sinusoidal bed?
% It is the depth proportional to the shallowest depth (riffle)

