% 2023-10-26 13:33:27.929355172 +0200
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
syms a1 a2 x z(x) f(x)
%ode = (0 == z + a1*diff(z,x) + a2*diff(z,x,2)) 
ode = a2*diff(z,x,2) + a1*diff(z,x) + z(x)==0;
zs = dsolve(ode)
%,z(0)==1);

% solution depends on if a1^2 - 4*a2 > 0
% (a1^2 - 4*a2)^(1/2)))/(2*a2)
r = a1/(2*a2);
%ri = 1i*sqrt(a1^2 - 4*a2)/(2*a2)
o = sqrt(4*a2 - a1^2)/(2*a2)
%rfun=matlabFunction(r); rfun(0.1,0.1)

if (0)
%C1*exp(-(x*(rr - ri)) - exp(-(x*(rr + ri))*(C1 - 1)
%C1*exp(-(x*(rr - ri)) - exp(-(x*(rr + ri))*C1 - exp(-(x*(rr + ri) 
%C1*exp(-x*(rr))*(exp(x*ri) - exp(-x*ri))  - exp(-(x*(rr + ri)
%C1*exp(-x*(rr))*(exp(x*ri) - exp(-x*ri))  - exp(-(x*(rr + ri)

% ReC1 = ReC2
% ImC1 = -ImC2
syms C1 C2
C1*exp(-(x*(a1 + (a1^2 - 4*a2)^(1/2)))/(2*a2)) + C2*exp(-(x*(a1 - (a1^2 - 4*a2)^(1/2)))/(2*a2))
C1*(exp(-(x*(a1 + (a1^2 - 4*a2)^(1/2)))/(2*a2)) + exp(-(x*(a1 - (a1^2 - 4*a2)^(1/2)))/(2*a2)))
C1*exp(-x*r)*(exp(-x*o*1i) + exp(+x*o*1i))
C1*exp(-x*r)*2*cos(x*o);
C1*exp(-x*r)*(exp(-x*o*1i) - exp(+x*o*1i))
C1*exp(-x*r)*(-2i)*sin(x*o);
%(exp(-x*ri) + exp(+x*ri))
end

% syms r o
syms c1 c2
zs = exp(-x*r)*(c1*cos(x*o) + c2*sin(x*o))


es = a2*diff(zs,x,2) + a1*diff(zs,x,1) + zs;

e0 = subs(es,x,0)

%syms b2 b4 z(x)
%ode2 = b4*diff(z,x,4) + b2*diff(z,x,2) + z == 0
%tf = dsolve(ode2,z)

% continuous transform:
a2*diff(z,x,2) + a1*diff(z,x) + z(x) == f(x);
(a2*(i o)^2 + a1*(i o) + 1) z == x(o)
(-a2*o^2 + a1*1i*o + 1) z == x(o)
response = z/x = 1/(1 + a1 i o - a2*o^2)

% impulse response
syms a1 a2 x y(x); dsolve(y + a1*diff(y,x) + a2*diff(y,x,2) == dirac(x))

