% Thu 23 Mar 13:59:46 CET 2017

syms a0 a1 a2 b1 b2 x dp

A = rand(5,1)
for idx=2
	f = (a0 + a1*sin(x) + a2*sin(2*x+dp))^idx
	%f = (a0 + a1*sin(x) + b1*cos(x) + a2*sin(2*x) + a3*cos(2*x))^idx
	%simplify(f2)
	double(subs(subs(subs(subs(subs(f,a0,A(1)),a1,A(2)),a2,A(3)),dp,A(4)),x,A(5)))
	f = expand(f);
	f = combine(f,'sincos')
	double(subs(subs(subs(subs(subs(f,a0,A(1)),a1,A(2)),a2,A(3)),dp,A(4)),x,A(5)))
end
children(f).'
% 
% f2 = a0^2 - (a2^2*cos(4*x))/2 - (a1^2*cos(2*x))/2 + a1^2/2 + a2^2/2 + 2*a0*a1*sin(x) - a1*a2*cos(3*x) + 2*a0*a2*sin(2*x) + a1*a2*cos(x)
% f3 = (3*a1^3*sin(x))/4 - (a1^3*sin(3*x))/4 + (3*a2^3*sin(2*x))/4 - (a2^3*sin(6*x))/4 + (3*a0*a1^2)/2 + (3*a0*a2^2)/2 + a0^3 + 3*a0^2*a1*sin(x) + (3*a1*a2^2*sin(x))/2 - (3*a0*a1^2*cos(2*x))/2 - (3*a0*a2^2*cos(4*x))/2 + 3*a0^2*a2*sin(2*x) + (3*a1^2*a2*sin(2*x))/2 + (3*a1*a2^2*sin(3*x))/4 - (3*a1^2*a2*sin(4*x))/4 - (3*a1*a2^2*sin(5*x))/4 + 3*a0*a1*a2*cos(x) - 3*a0*a1*a2*cos(3*x)


