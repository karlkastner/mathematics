
% dy/dt = f(y) + g(y)

%-> simpler (see randall leveque)
%e(A+B) ~ I + (A+B)/2 + 1/4*(A+B)^2
%e(A+B) ~ e(A/2)*e(B)*e(A/2)
%       ~ (I+A/2+A^2/4)*(I + B + B^2/2)*(I + A/2 + A^2/4)
%       ~ (I + A/2 + A^2/4 + B + 1/2*A*B + 1/4*A*B^2 + 

syms f0 g0 dfdy dgdy d3ydt3 d2ydt2

% negect second derivatives here
f = @(y) f0 + dfdy*(y-y0)
g = @(y) g0 + dgdy*(y-y0)

dydt = f0 + g0

% with forward euler as substages:
% this is not second order accurate bc the substage scheme is only first order:
ya = y0 + 1/2*dt*f(y0)
yb = ya + dt*g(ya)
yc = yb + 1/2*dt*f(yb)

% with heuns second order method for substages:
yap = y0 + 1/2*dt*f(y0)
ya  = y0 + 1/4*dt*(f(y0) + f(yap)) %y0 + 1/2*dt*f(y0)))
ybp = ya + dt*g(ya)
yb  = ya + dt*1/2*(g(ya) + g(ybp))
ycp = yb + 1/2*dt*f(yb)
yc  = yb + dt*1/4*(f(yb) + f(ycp))

%ynext = y0 + dt*dy0dt + 1/2*dt^2*d2y0dt2
%ynext = y0 + dt*(f0 + g0) + 1/2*dt^2*(df/dt + dg/dt)
ynext = y0 + dt*(f0 + g0) + 1/2*dt^2*(dfdy*dydt + dgdy*dydt) + 1/6*d3ydt3
%ynext = y0 + dt*(f0 + g0) + 1/2*dt^2*(d2ydt2) + 1/6*d3ydt3
e=expand(simplify(yc-ynext))
eleading = simplify(subs(subs(subs(e,dt^6,0),dt^5,0),dt^4,0))

if (0)
% example that without two half steps it is only first order accurate:
yap = y0 + dt*f(y0)
ya  = y0 + 1/2*dt*(f(y0) + f(yap)) %y0 + 1/2*dt*f(y0)))
ybp = ya + dt*g(ya)
yb  = ya + dt*1/2*(g(ya) + g(ybp))
e2=expand(simplify(yb-ynext))
eleading2 = simplify(subs(subs(subs(e2,dt^6,0),dt^5,0),dt^4,0))
end

% for error estimate two half vs full step
yap = y0 + 1/2*dt*f(y0)
ya  = y0 + 1/4*dt*(f(y0) + f(yap)) %y0 + 1/2*dt*f(y0)))

yaap = y0 + 1/4*dt*f(y0)
yaa  = y0 + 1/8*dt*(f(y0) + f(yaap))
yabp = yaa + 1/4*dt*f(yaa)
yab  = yaa + 1/8*dt*(f(yaa) + f(yabp))

e = ya-yab;
expand(simplify(e))
