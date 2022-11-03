% Thu 17 Nov 15:18:24 CET 2016
pflag = 0;

n  = 500;
nf = round(n/3);

x0 = (0:n/4-1)'/(n);
y0 = 1-8*x0;
x = (0:n/2-1)'/(n); %(linspace(0,1,n)';
y1 = -cos(2*pi*x);
y2 = -sign(y1);
y = [y0;y1;y2];

n = length(y);
x = (0:n-1)'/(n);

mu = 0; sd = 0.0;
sk=0;k=100000;
r = 0; %pearsrnd(mu,sd,sk,k,n,1);
r = sd*randn(size(y));
y = flipud(y);
y0 = y;
y = y + r;

wh  = hanwin((1:nf)');
wk  = kaiserwin((1:nf)');
wl = lanczoswin((1:nf)');
wr = rectwin((1:nf)');
wt = triwin((1:nf)');
	
%w = ones(nf,1);

p=0.16;

figure(1);
clf

ll = 0;
if (pflag) figure(1); clf; subplot(3,1,1); else subplot(3,1,1); end

d = 0;

ll = 1;
d = 0;
y_ = ll*wmedfilt(wh,y,d);
plot(x,y,'.');
hold on
plot(x,y0,'k');
plot(x,[  wordfilt(wh,y-y_,1-p,d)+y_, ...
	  ... %wordfilt(w,y-y_,0.5,d)+y_, ...
	  wmedfilt(wh,y,d), ...
	  wmedfilt(wl,y,d), ...
	  wordfilt(wh,y-y_,p,d)+y_, ...
	  wmeanfilt(wh,y,d) ...
	  wmeanfilt(wl,y,d) ...
		]) 
ylabel('y')
xlabel('x')
legend('location','eastoutside','samples','original','p16', 'me han', 'me lanczos', 'p84','mean')
if (pflag)
	pdfprint(1,'/home/pia/phd/src/test/img/order-filter-a.pdf',1,1/1,'svg');
end

if (pflag) figure(2); clf; subplot(3,1,1); else subplot(3,1,2); end
ll=0;
y_ = ll*wmedfilt(wh,y,d);
plot(x,y,'.');
hold on
plot(x,y0,'k');
plot(x,[  wordfilt(wh,y-y_,1-p,d)+y_, ...
	  ... %wordfilt(w,y-y_,0.5,d)+y_, ...
	  wmedfilt(wh,y,d), ...
	  wordfilt(wh,y-y_,p,d)+y_, ...
	  wmeanfilt(wh,y,d) ...
		]) 
ylabel('y')
xlabel('x')
if (pflag)
	pdfprint(2,'/home/pia/phd/src/test/img/order-filter-b.pdf',1,1/1,'svg');
end




if (pflag) figure(3); clf; subplot(3,1,1); else subplot(3,1,3); end
d = 1;
ll=1;
y_ = ll*wmedfilt(wh,y,d);
plot(x,y,'.');
hold on
plot(x,y0,'k');
plot(x,[  wordfilt(wh,y-y_,1-p,d)+y_, ...
	  ... %wordfilt(w,y-y_,0.5,d)+y_, ...
	  wmedfilt(wh,y,d), ...
	  wmedfilt(wl,y,d), ...
	  wmedfilt(wr,y,d), ...
	  wmedfilt(wt,y,d), ...
	  wordfilt(wh,y-y_,p,d)+y_, ...
	  wmeanfilt(wh,y,d) ...
		]) 
legend('sample','orig','d16','me','l','r','t','d84','mean')
xlabel('x')
ylabel('y')
if (pflag)
	pdfprint(3,'/home/pia/phd/src/test/img/order-filter-c.pdf',1,1/1,'svg');
	pdfprint(1,'/home/pia/phd/src/test/img/order-filter.pdf',1,1/1,'svg');
end
