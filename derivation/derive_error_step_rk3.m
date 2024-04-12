

syms y0 d1y d2y d3y d4y y5 dt e

ymid  = y0 + 0.5*dt*d1y
d1ymid = d1y + 0.5*dt*d2y
y1    = y0 - dt*d1y + 2*dt*d1ymid
d1y1    = d1y - dt*d2y + 2*dt*(d2y + 0.5*dt*d3y)
ynext_num = y0 + dt*(d1y + 4*d1ymid + d1y1)/6

ynext  = y0  + d1y*dt + 1/2*d2y*dt^2 + 1/6*d3y*dt^3 + 1/24*d4y*dt^4 % + 1/120*y5*dt^5

%dynext = d1y + d2y*dt + 1/2*d3y*dt^2 + 1/6*d4y*dt^3 + 1/24*y5*dt^4
eq = ynext == ynext_num + e
solve(eq,e)

