% Wed Jun 13 20:29:19 MSK 2012
% Karl Kästner, Berlin

syms x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4 x y z

A = [	1 1 1 1;
	x1 x2 x3 x4;
	y1 y2 y3 y3;
        z1 z2 z3 z4 ]

p = [1; x; y; z]

d1 = [x1 - x; y1 - y; z1 - z]
d2 = [x2 - x; y2 - y; z2 - z]
d3 = [x3 - x; y3 - y; z3 - z]
d4 = [x4 - x; y4 - y; z4 - z]
d1 = d1.'*d1
d2 = d2.'*d2
d3 = d3.'*d3
d4 = d4.'*d4

% find test function coefficients
C = inv(A) %*det(A)

c1 = (C(1,:)*p)
c2 = (C(2,:)*p)
c3 = (C(3,:)*p)
c4 = (C(4,:)*p)

% maximise error
f = simple(d1*c1 + c2*d2 + c3*d3 + c4*d4)
dx = diff(f,x)
dy = diff(f,y)
dz = diff(f,z)
sx = simple(solve(dx,x))
sy = simple(solve(dy,y))
sz = solve(dy,z)
sz = simple(sz)
f_ = subs(subs(subs(f,x,sx),y,sy),z,sz)
simple(f_)

