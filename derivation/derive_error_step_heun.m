

syms y y1 y2 y3 y4 y5 dt e

y_p_num   = y   + dt*y1 
dy_p_num  = y1  + dt*y2 
%dy_pre = y1 + y2*dt + 1/2*y3*dt^3 + 0*1/6*y3*dt^4
ynext_num  = y + 1/2*dt*(y1 + dy_p_num)
ynext  = y  + y1*dt + 1/2*y2*dt^2 + 1/6*y3*dt^3 + 1/24*y4*dt^4 % + 1/120*y5*dt^5

%dynext = y1 + y2*dt + 1/2*y3*dt^2 + 1/6*y4*dt^3 + 1/24*y5*dt^4
eq = ynext == ynext_num + e
solve(eq,e)

