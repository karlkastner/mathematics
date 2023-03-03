% Thu 22 Mar 11:19:13 CET 2018

n = 1e6;


%dist = 'logu';
dist = 'logtri';


switch (dist)
case {'logtri'}
a = 2;
%b = 3.5;
%b = a*(1+sqrt(eps));
c = 5;
b = a*(1+sqrt(eps));
%b = c*(1-sqrt(eps));

%	funrnd = @(n) zeros(n,1);
%	funcdf = @(x) zeros(size(x));
	funinv = @(F) logtriinv(a,b,c,F);
	funrnd = @(n) logtrirnd(a,b,c,n,1);
	funcdf = @(x) logtricdf(a,b,c,x);
	funpdf = @(x) logtripdf(a,b,c,x);
case {'logu'}
a = 2;
%b = 4;
c = 5;
b = c-sqrt(eps);
	funrnd = @(n) trirnd(a,b,c,n,1);
	funcdf = @(x) tricdf(a,b,c,x);
	funpdf = @(x) tripdf(a,b,c,x);
%else
	a = 2;
	b = 5;
	funrnd = @(n) logurnd(a,b,n,1);
	funpdf = @(x) logupdf(a,b,x);
	funcdf = @(x) logucdf(a,b,x);
end


r = funrnd(n);
r = real(r);
'mean'
mean(r)
meanlogu(a,b)
'var'
var(r)
varlogu(a,b)

pfun = @semilogx;
%pfun = @plot

figure(1);
clf();

subplot(2,2,1)
% = linspace(0,10,19)';
x = linspace(1,10,101)';
%[h,x] = hist(r,x);
h = histc(r,x);
h = h/(n*(x(2)-x(1)));
h(:,2) = funpdf(x);
%x = centre(x); 
%bar(x,h)
pfun(x,h);
title('pdf');
hold on
vline([a b c])

disp('int f == 1')
sum(h)*(x(2)-x(1))
%pause

subplot(2,2,2);
H = cumintL(h(:,1),x);
H(:,2) = funcdf(x);
pfun(x,H);
%plot(x,cumintL(h,x));
title('cdf');
vline([a b c])

subplot(2,2,3)
% = linspace(0,10,19)';
m = 101;
x_ = logspace(log10(a)-0.05,log10(c)+0.05,m)';
%[h,x] = hist(r,x);
h = histc(r,x_);
h = h*m/n;
%x = centre(x); 
bar(x_,h)
set(gca,'xscale','log')
xlim(limits(x_));
title('Exponentially spaced bins');

subplot(2,2,4)
cla
n = 100;
F_ = (1:n)'/(n+1);
x_ = funinv(F_);
x_ = real(x_);
plot(H(:,2),x);
hold on
plot(F_,x_,'.');
hline([a b c])

