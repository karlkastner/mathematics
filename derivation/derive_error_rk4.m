% Tue 12 Mar 15:46:41 CET 2024
syms y  y0 aa h
ydot = @(t,y) aa*y

t = 0;
yh = step_rk4(ydot,y0,t,h) 
yh = expand(yh)
% y(t)  = y0*exp(a*t)
%yh_ = y0*(1 + a*h + 1/2*a^2*h^2)

yh_ = taylor(y0*exp(aa*h),h,0,'order',6)
%(aa^5*y^5)/120 + (aa^4*y^4)/24 + (aa^3*y^3)/6 + (aa^2*y^2)/2 + aa*y + 1
 
e = yh - yh_
% so the truncation error is equal to -h^5*d^5y/dt

