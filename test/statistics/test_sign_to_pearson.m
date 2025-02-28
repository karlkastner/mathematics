% 2024-12-17 23:50:44.379783827 +0100
% Karl Kastner, Berlin

r = 0.9;
 n = 1e6;
 x = randn(n,1);
 y = r*x + sqrt(1-r^2)*randn(n,1);
 c=corr(sign(x),sign(y))
 r_ = sign_to_pearson(r)
 c_ = pearson_to_sign(c)
% c_from_r = 2/pi*asin(r)
% sin(pi/2*c)
% 2/pi*asin(c)
% pi/2*sin(r)
