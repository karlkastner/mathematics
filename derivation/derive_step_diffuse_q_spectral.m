F(I - e*D2)
 (F - e*F*D2)

D2 = [1,-2,1]
Ansatz: V*E*V'
% set x = 0
[sin(o*(x-1/n)) + cos(o*(x-1/n)) + -2*(sin(o*x) + cos(o*x)) + sin(o*(x+1/n) + cos(o*(x+1/n))] = e*[sin(o*x) + cos(o*x)]
[sin(o*(-1/n)) + cos(o*(-1/n)) + -2*1 + sin(o*(1/n) + cos(o*(1/n))] = e*1
[+ cos(o*(-1/n)) + -2*1 + cos(o*(1/n))] = e
e = 2*(-1 + cos(o*1/n))
o = 2*pi*(0:n/2);
e = 2*(cos(o)-1)

[a*sin(o*(x-1/n)) + b*cos(o*(x-1/n)) + -2*(a*sin(o*x) + b*cos(o*x)) + a*sin(o*(x+1/n) + b*cos(o*(x+1/n))] = e*[a*sin(o*x) + b*cos(o*x)]
% x = pi/2
[a*cos(o*(-1/n)) - b*sin(o*(-1/n)) + -2*a + a*cos(o*(1/n) - b*sin(o*(1/n))] = e*[a*cos(o*x)]
[a*cos(o*(-1/n)) -2*a + a*cos(o*(1/n))] = e*[a*cos(o*x)]

% general
[a*(sin(o*(x-1/n)+sin(o*(x+1/n))) + b*(cos(o*(x-1/n))+cos(o*(x+1/n))) + -2*(a*sin(o*x) + b*cos(o*x))] = e*[a*sin(o*x) + b*cos(o*x)]
[2*b*cos(o*x)*cos(o*1/n) + 2*a*sin(o*x)*cos(o*1/n) + -2*a*sin(o*x) - 2*b*cos(o*x)] = e*[a*sin(o*x) + b*cos(o*x)]

compare coeffs
[2*b*cos(o*x)*cos(o*1/n) -2*b*cos(o*x)] = e*[b*cos(o*x)]
[2*a*sin(o*x)*cos(o*1/n) -2*a*sin(o*x)] = e*[a*sin(o*x)]
Duplicate eigs, except 0
% note that eigenvectors are not unique bc of duplicate eigenvalues!

e = 2*(cos(o/n)-1)
[-2*a*sin(o*x)] = e*[a*sin(o*x)]
2*b*cos(o*1/n) + b = e*b
-2*a = e*a -> a = 0
o = 2*pi*[(0:n/2),(n/2-1)-1:1]; e = 2*(cos(o/n)-1);
e_ = diag(v'*D2*v)' 
v=[cos(x*(0:n/2))-1i*sin(x*(0:n/2))]; v=[v,fliplr(conj(v(:,2:end-1)))];



(I - e*D2)^-1*x
F^-1 F (I - e*D2)^-1 F^-1*F*x)
F^-1 F ( F V( I - e*E)*V' )^-1 *F*x)

V'*(I - e*V*E*V')*V

V'*dz/dt = V'*(I - e*V*E*V')*z
dV'z/dt = V'*(I - e*V*E*V')*V*V'z
d tz/dt = (I - e*E)* tz
tz0 = V*z0
tz = exp(-(1-e*E))*tz0
z = V*tz

B*F = V'
B   = V'*(F)^-1
F*C = V
C = n*B'

tilde z = V'z
dtilde z/dt = (I - e*E) tilde z


[ax*sin(o*(x-d))+ay*sin(o*(y-d)) + bx*cos(o*(x-d))+by*cos(o*y-d))
 -4*(ax*sin(o*(x))+ay*sin(o*(y)) + bx*cos(o*(x))+by*cos(o*y))
  ax*sin(o*(x-d))+ay*sin(o*(y-d)) + bx*cos(o*(x-d))+by*cos(o*y-d))] = e*(ax*sin(ox) + ay*sin(oy) + bx*cos(ox) + by*cos(oy))
set x = 0, y=0

% bx=by = 1
[bx*ex*cos(ox*d)+by*ey*cos(oy*d))
 -2*(ex*bx+ey*by)
  bx*ex*cos(ox*(d))+by*ey*cos(oy*d))] = e*(bx + by)
  2*bx*ex*cos(ox*d)+2*by*ey*cos(oy*d) - -2*(ex*bx+ey*by) = 2*(bx + by)

% first derivative
[-a*sin(o x-d) - b*cos(o x-d) + a*sin(o*x) + b*cos(o*x)] = e*(a*sin(ox) + b*cos(ox))
% x = 0 
% set a = b = 1
[-a*sin(-d) - b*cos(-d) + b] = e*(b)
[-1i*sin(-d) - cos(-d) + 1] = e

% other direction
[a*sin(o x+d) + b*cos(o x+d) - a*sin(o*x) - b*cos(o*(x))] = e*(a*sin(ox) + b*cos(ox))
% x =0
[a*1i sin(o d) + b*cos(d) -  b] = e*b
[1i sin(o d) + cos(o*d) -  1] = e*b



