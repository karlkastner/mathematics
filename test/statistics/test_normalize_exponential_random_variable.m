% Fri 24 Mar 15:48:04 CET 2023

x = exprnd(1,1e7,1);
 x=[x,normalize_exponential(x)];
 [h,xi]=hist((x-mean(x))./std(x),innerspace(-3,3,20));
 h=h./(sum(h)*(xi(2)-xi(1)));
 clf;
 bar(xi,h);
 hold on;
 plot(xi,normpdf(xi));
 skewness(x), kurtosis(x)

