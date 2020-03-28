% 2012 Apr 26 15:26 MCT
% Karl Kästner, Berlin

syms x1 x2 x3 y1 y2 y3 x y

A = [	1 1 1;
	x1 x2 x3;
	y1 y2 y3 ]

p = [1; x; y]

d1 = [x1 - x; y1 - y]
d2 = [x2 - x; y2 - y]
d3 = [x3 - x; y3 - y]
d1 = d1.'*d1
d2 = d2.'*d2
d3 = d3.'*d3

% find test function coefficients
C = inv(A) %*det(A)

c1 = (C(1,:)*p)
c2 = (C(2,:)*p)
c3 = (C(3,:)*p)

% maximise error
f = simple(d1*c1 + c2*d2 + c3*d3)
dx = diff(f,x)
dy = diff(f,y)
sx = solve(dx)
sy = solve(dy,y)
f_ = (subs(subs(f,x,sx),y,sy))
simple(f_)

%/det(A))
%- ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)))

