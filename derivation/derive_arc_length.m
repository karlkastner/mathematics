% 2016-04-14 17:28:38.020451648 +0200
% Karl Kastner, Berlin

%sys a b t
syms x1 x2 x3 x4 y1 y2 y3 y4 a1 a2 a3 a4 t


% this does not seem to be more accurate, because t is not linear
A = [1 -1  1; 
     1  0  0
     1  1  1];
x = [x1; x2; x3];
y = [y1; y2; y3;];
cx = A \ x
cy = A \ y

f  = (cx(2) + cx(3)*t)^2 + (cy(2) + cy(3)*t)^2
f  = a + b*t + c*t^2
f  = collect(f,t)
I  = int(sqrt(f),t),
I1 = (subs(I,t,0)-subs(I,t,-1))
I2 = (subs(I,t,1)-subs(I,t,0))

%I1 = 1/b*( (a*abs(a)) - (a - b/2)*abs(a - b/2) )
%I2 = 1/b*( (a + b/2)*abs(a + b/2) - (a*abs(a)) )

%I1 = (a^(1/2)*b)/(4*c) - (b/(4*c) - 1/2)*(a - b + c)^(1/2) + (log(b/(2*c^(1/2)) + a^(1/2))*(- b^2/4 + a*c))/(2*c^(3/2)) - (log((b/2 - c)/c^(1/2) + (a - b + c)^(1/2))*(- b^2/4 + a*c))/(2*c^(3/2))
%I2 = (b/(4*c) + 1/2)*(a + b + c)^(1/2) + (log((b/2 + c)/c^(1/2) + (a + b + c)^(1/2))*(- b^2/4 + a*c))/(2*c^(3/2)) - (a^(1/2)*b)/(4*c) - (log(b/(2*c^(1/2)) + a^(1/2))*(- b^2/4 + a*c))/(2*c^(3/2))

if (0)
A = [1 -1 1 -1;
     1 -1/2 1/4 -1/8;
     1 1/2 1/4 1/8;
     1  1 1 1];

x = [x1; x2; x3; x4];
y = [y1; y2; y3; y4];
cx = A \ x
cy = A \ y

%f  = (cx(2) + cx(3)*t)^2 + (cy(2) + cy(3)*t)^2 
f  = a1 + a2*t + a3*t^2 + a4*t^3;
f  = collect(f,t)
I  = int(sqrt(f),t),
I1 = (subs(I,t,0)-subs(I,t,-1))
I2 = (subs(I,t,1)-subs(I,t,0))

end

