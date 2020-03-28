% 2012-07-05 17:14:28 UTC
% Karl Kästner, Berlin

syms h x y
f = exp(-sqrt(x^2 + y^2))      
df1 = diff(f,x)
df2 = diff(df1,x)
df3 = diff(df2,x)
df4 = diff(df3,x)

df1 = simple(df1/f)
df2 = simple(df2/f)
df3 = simple(df3/f)
df4 = simple(df4/f)

df1
df2
df3
df4

%df4_ = (-12*x^2*y^2 + 3*y^4  + x^6 + x^4*y^2)/(x^2 + y^2)^3 + (-5*x^6 - 11*x^4*y^2 - 12*x^2*y^2 - 6*x^2*y^4 + 3*y^4)/(x^2 + y^2)^(7/2)
%(-6*x^2 + 3)*(x^2 + y^2)^2 + (-18*x^2 + 6*x^4 )*(x^2 + y^2) + 15*x^4
%3*(x^2 + y^2)^2 + (-18*x^2 + x^4)*(x^2 + y^2) + 15*x^4
%df4_ = -3*y^2*(2*x^4 + 2*x^2*y^2 + 4*x^2 - y^2)/(x^2 + y^2)^(7/2) + (x^6 + x^4*y^2 - 12*x^2*y^2 + 3*y^4)/(x^2 +y^2)^3
df1 = subs(subs(df1,x,h),y,h)
df2 = subs(subs(df2,x,h),y,h)
df3 = subs(subs(df3,x,h),y,h)
df4 = subs(subs(df4,x,h),y,h)


df1_ = -1/(2^(1/2))
df2_ = 1/2 - 1/(2^(3/2)*h)
df3_ = 3/(4*h) + ((3*h^-2 - 2))/(2^(5/2))
df4_ = (- 12*h^-1 - 9*h^-3 - 9*2^(1/2)*h^-2 + 2*2^(1/2) )/(8*2^(1/2))

subs(df1_ - df1,h,0.111)
subs(df2_ - df2,h,0.111)
subs(df3_ - df3,h,0.111)
subs(df4_ - df4,h,0.111)


df1_ = -1/(2^(1/2))
df2_ = -1/(2*2^(1/2)*h)
df3_ = ((3*h^-2))/(4*2^(1/2))
df4_ = (- 9*h^-3 )/(8*2^(1/2))

