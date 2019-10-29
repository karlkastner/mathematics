syms x1 x2 x3
dx1 = x1 - x1*x2 - x2^3 + x3*(x1^2 + x2^2 - 1 - x1 + x1*x2 + x2^3)
dx2 = x1 - x3*(x1 - x2 + 2*x1*x2)
dx3 = (x3 - 1)*(x3 + 2*x3*x2^2 + x3^3)

s3=solve(dx3,x3)

solve(subs(dx2,x3,s3(1)),x2)
s2 = solve(subs(dx2,x3,s3(2)),x2)
solve(subs(dx2,x3,s3(3)),x2)
solve(subs(dx2,x3,s3(4)),x2)

solve(subs(subs(dx1,x3,s3(2)),x2,s2),x1)
subs(subs(dx1,x3,s3(2)),x2,s2)

d = dx1 + dx2 + dx3
s = solve(d,x1)

d_ = subs(d,s(1))
%olve(d_,d_)

subs(dx1,x3,1)
subs(dx2,x3,1)

subs(subs(subs(dx1,x3,1),x2,sqrt(3/4)),x1,0.5)
subs(subs(subs(dx2,x3,1),x2,sqrt(3/4)),x1,0.5)
subs(subs(subs(dx3,x3,1),x2,sqrt(3/4)),x1,0.5)
