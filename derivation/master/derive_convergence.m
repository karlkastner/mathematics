% 2012-05-30 09:48:55 UTC
% Karl Kästner, Berlin

% this is in support of the FDM chapter of the thesis,
% which attempts to show, that the finite difference scheme (slowly) converges
% at the origin of the Coulomb potential
% the derivation (but not the finding) was challenged at the defence by Prof. 

syms x y h

f1 = -x/(exp((x^2 + y^2)^(1/2))*(x^2 + y^2)^(1/2))
 
f2 = (x^2/(x^2 + y^2) - 1/(x^2 + y^2)^(1/2) + x^2/(x^2 + y^2)^(3/2))/exp((x^2 + y^2)^(1/2))
 
f3 = ((3*x - x^3)/(x^2 + y^2)^(3/2) + (3*x)/(x^2 + y^2) - (3*x^3)/(x^2 + y^2)^2 - (3*x^3)/(x^2 + y^2)^(5/2))/exp((x^2 + y^2)^(1/2))
 
f4 = -((6*x^2 - 3)/(x^2 + y^2)^(3/2) - 3/(x^2 + y^2) + (18*x^2 - x^4)/(x^2 + y^2)^2 + (18*x^2 - 6*x^4)/(x^2 + y^2)^(5/2) - (15*x^4)/(x^2 + y^2)^3 - (15*x^4)/(x^2 + y^2)^(7/2))/exp((x^2 + y^2)^(1/2))

f1 = simple(subs(subs(f1/f0,x,h),y,h))
f2 = simple(subs(subs(f2/f0,x,h),y,h))
f3 = simple(subs(subs(f3/f0,x,h),y,h))
f4 = simple(subs(subs(f4/f0,x,h),y,h))

f1
f2
f3
f4

