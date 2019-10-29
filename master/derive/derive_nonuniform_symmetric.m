% 2012-01-21 21:26:56 UTC
% Karl Kästner, Berlin

syms fp1 f fm1
%syms h hm1 hp1
syms xm1 x xp1

F = [fm1 f fp1].'
hm1 = x - xm1
h   = xp1 - x

eq1 = -1/2*(2*(fp1 -f)/(h*(h+hm1)) - 2*(f - fm1)/(hm1*(h+hm1)))
%eq2 = fm1 * -1/(hp1*(h+hm1)) + fp1 * -1/(hp1*(h+hm1)) - f * (-1/(hp1*(h+hm1)) + -1/(hp1*(h+hm1)))
Ls = (h + hm1)/2
eq3 = 1/Ls*[ -1/(2*hm1), (h+hm1)/(2*h*hm1), -1/(2*h)]
eq3_s = [ -1/(2*hm1), (h+hm1)/(2*h*hm1), -1/(2*h)]

%expand(eq1 - eq2)
%expand(eq1 - eq3)
%expand(eq2 - eq3)
[n1 d1] = numden(eq1);
[n3 d3] = numden(eq3)
[n3 d3] = numden(eq3*F)
n1 - n3
d1 - d3
[n3_ d3_] = numden(eq3_s)
[n3_ d3_] = numden(eq3_s*F)
%%%
syms xm2 xm1 x0
A = [ (x1 - x0) (x1  - x0)^2/2;
      (xm1 - x0) (xm1 - x0)^2/2] 
B = [ 0 -1 1;
      1 -1 0]
s = A \ B
s_2 = s(2,:)

syms xm2 xm1 x0
syms hll hl hr hrr 
A = [ (x2 - x0)  (x2  - x0)^2/2  (x2 - x0)^3/6 (x2 - x0)^4/24
      (x1 - x0)  (x1  - x0)^2/2  (x1 - x0)^3/6 (x1 - x0)^4/24;
      (xm1 - x0) (xm1 - x0)^2/2  (xm1 - x0)^3/6 (xm1 - x0)^4/24;
      (xm2 - x0) (xm2 - x0)^2/2  (xm2 - x0)^3/6 (xm2 - x0)^4/24]
% = [ hrr hrr^2/2 hrr^3/6 hrr^4/24;
%     hr  hr^2/2  hr^3/6  hr^4/24;
%     -hl  hl^2/2  -hl^3/6  hl^4/24;
%     -hll hll^2/2 -hll^3/6 hll^4/24]
 
B = [ 	0 0  -1 0 1;
	0 0  -1 1 0;
        0 1  -1 0 0;
        1 0  -1 0 0]
s = A \ B
s_4 = s(2,:)

err =simple(s_4 - [0 0 s_2 0 0])

