% Sat 10 Dec 21:05:15 CET 2016
% Karl Kastner, Berlin
%
%% arc length of a two dimensional curve
%%
%% 8th order accurate
%% does not require the segments length to vary smoothly
%%
%% note: the curve can be considered parametric, e.g. x = x(t), y=y(t) and
%% 	and t = t(s), but the error term contains derivatives of t,
%%      thus a non smooth t (strongly varying distance between points)
%%      requires the scaling as done below
%
function S = arclength(x,y,order)
	if (nargin() < 3)
		order = 5;
	end

	% 1st approximate segment length by cord length
	x  = rvec(x);
	y  = rvec(y);
	dx = diff(x);
	dy = diff(y);
	% approximate arc length
	ds = hypot(dx,dy);
	% chord length, second order accurage
	S  = [0, cumsum(ds)];

if (order > 2)
	% it is important not to use point indices to set-up the polynomial,
	% but to scale the segments with the chord length
	% only this gives higher order accurcary when dS is not constant

	% spline (picewise hermite polynomial)
	px  = spline(S,x);
	py  = spline(S,y);

	% derivative polynomial
	dpx = fnder(px)
	dpy = fnder(py);

	% evaluate derivative at points and mid-points
	% interpolate curve to mid sections
	Si(1:2:2*length(S)-1) = S;
	Si(2:2:end-1)         = 0.5*(S(1:end-1)+S(2:end));
	dx  = ppval(dpx,Si);
	dy  = ppval(dpy,Si);

	% arc length function
	f = hypot(dx,dy);

	% evaluate with simpsons's rule
	S = [0,1/6*cumsum( ds.*(f(1:2:end-2)+4*f(2:2:end-1)+f(3:2:end)) )];
	%F = [0,1/3*cumsum( ds.*(f(1:2:end-2)+4*f(2:2:end-1)+f(3:2:end)) )];
end

if (0)

	% S = (0:n-1)/(n-1);
	% Si = (0:2*(n-1))/(2*(n-1));
	x = interp1(S,x,Si,'spline');
	y = interp1(S,y,Si,'spline');

	% aprroximate curve in each section as a quadratic polynomials in x and y
	A = [1 -1 1;
             1  0 0;
             1  1 1];
	% coefficients
	cx = A \ [x(1:2:end-2);x(2:2:end-1);x(3:2:end)];
	cy = A \ [y(1:2:end-2);y(2:2:end-1);y(3:2:end)];
	% coefficients of derivative
	dcx = [cx(2,:); 2.*cx(3,:)];
	dcy = [cy(2,:); 2.*cy(3,:)];
	% value of the derivative at mid and end-points
	dx = A(:,1:2)*dcx;
	dy = A(:,1:2)*dcy;
dx
dy
pause
	% arc length function f = sqrt(dx^2+dy^2)
	f = hypot(dx,dy);
	% midpoint
%	F = [0,cumsum(f(2,:))];
	% trapezoidal rule
%	F = [0,0.5*cumsum(f(1,:)+f(3,:))];
	% simpsons rule
	F = [0,1/3*cumsum(f(1,:)+4*f(2,:)+f(3,:))];
end
end % arclength

