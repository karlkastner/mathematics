% 2016-02-29 13:36:34.401870431 +0100
% Karl Kastner, Berlin

n=100;
 p=(1:n)/(n+1);
 c=normpdf(norminv(p));
 q=(1-p);
  sd = 1./c.*sqrt(p.*q);
 plot(p,sd)

