% 2016-02-26 14:47:36.689926277 +0100

n = 100;
m = 1e4;
x = randn(n,m);

me = median(x);
p = 0.05;
q  = quantile(me,[p/2,0.5,1-p/2])
%q  = quantile(me,[p,0.5,1-p])
%p = 0.16;
%q  = quantile(me,[p,0.5,1-p])

disp('predicted:')
[qme, qs, ql, qr] = median_man(x,0.95);
disp('empirical:')
median(ql)
median(qr)


if (0)

[me s l u] = median_man((x),1-2*normcdf(-1));

a = [mean(l) quantile(me,normcdf([-1,1])), mean(u)]
a(1)./a(2)
a=[sqrt(mean(s.^2)) std(median(x)) 1.2533/sqrt(n),1.2533/sqrt(n)*(1-(4-pi)/(2*n))]
a(2)./a(1)

[me s l u] = median_man((x),1-2*normcdf(-2));
[mean(l) quantile(me,normcdf([-2,2])), mean(u)]

%serr(l)
%serr(u)

% cadwell correction of quantiles
% TODO in median ci: get pdf (c0) and cdf (p) of l and u, scale both by the value
n=1000; x = rand(n,1e5); mean(quantile(x,normcdf(-1)))./normcdf(-1), c0=normpdf(-1); p=normcdf(-1); 1/(n*c0^2/(n*c0^2 +p*(1-p)))
s=-1; n=10; x = rand(n,1e5); (mean(quantile(x,normcdf(s)))./normcdf(s)), c0=normpdf(-s); p=normcdf(-s); (n*c0^2/(n*c0^2 +p*(1-p)))

end

