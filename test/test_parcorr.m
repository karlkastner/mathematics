% Sat 25 Sep 12:29:50 CEST 2021
a=0.5;
b = 1-a;
%sqrt(1-a^2);
n = 1e5;
x=randn(n,1);
e   = 1e-3*randn(n,1);

if (0)

figure(1);
clf

n=5;
yi = x;
yis = x;
for idx=1:n
    yr = yi;
    yrs = yis;
    for jdx=1:n
	subplot(n,n,n*(idx-1)+jdx);
	parcorr(yr);
	title(['h' num2str(jdx-1), 'l', num2str(idx-1)]);
	yr = yr - filter(b,[1,-a],yr);
	if (mod(jdx,2) == 1)		
		yrs = yrs - fliplr(filter(b,[1,-a],fliplr(yrs)));
		yrs = yrs - filter(b,[1-a],yrs);
		p = parcorr(yr);
		hold on
		np = length(p);
		plot(0:np-1,p,'b*');
		axis auto
	end
    end
    yi = filter(sqrt(1-a^2),[1,-a],yi);
sys = arima(idx,0,0);
sys = estimate(sys,e(end/2:end)+yrs(end/2:end))
syms x;
vpa(expand((1-a*x^-1)^idx),2)
armax(yi(end/2:end)+e(end/2:end),[idx,0])
    if (mod(idx,2) == 1)		
		yis = filter(b,[1-a],fliplr(filter(b,[1,-a],fliplr(yis))));
sys = arima(2*idx,0,0);
sys = estimate(sys,e(end/2:end)+yis(end/2:end))
syms x;
vpa(expand(2*x^-1*(1-a*x^+1)*(1-a*x^-1)),2)
armax(yis(end/2:end)+e(end/2:end),[2*idx,1])
pause
    end
end

end

if (1)
figure(1);
clf()

if (0)
% one sided lowpass
y=(filter(sqrt(1-a^2),[1,-a],x));
subplot(2,4,1)
parcorr(y);

% one sided highpass
yhp=(x-filter(sqrt(1-a^2),[1,-a],x));
yhp_ = filter([1-sqrt(1-a^2),-a],[1,-a],x);
rms(yhp - yhp_)
subplot(2,4,2)
parcorr(yhp(:,1));

% one sided bandpass
ybp=filter(sqrt(1-a^2),[1,-a],(x-filter(sqrt(1-a^2),[1,-a],x)));
ybp_ = filter(sqrt(1-a^2)*[1-sqrt(1-a^2),-a],[1,-2*a,a^2],x);
rms(ybp-ybp_)
subplot(2,4,3)
%plot(y)
parcorr(ybp);

% second-order lowpass
subplot(2,4,4)
y2=filter(sqrt(1-a^2),[1,-a],filter(sqrt(1-a^2),[1,-a],x));
y2_=filter(1-a^2,[1,-2*a,a^2],x);
rms(y2-y2_)
%filter([1,-a],[1,-a],x);
%y(:,1) = yhp;
parcorr(y2)
end

a = 0.9;
b = sqrt(1-a^2);
b = 1-a;

x = zeros(length(x),1);
x(end/2) = 1;

if (1)
% two sided lowpass
y=lowpass1d_implicit(cvec(x),a);
y2 = fliplr(filter(b,[1,-a],fliplr(rvec(x))));
y2 = filter(b,[1,-a],rvec(y2));
y(:,2)=y2;
%filter(b,[1,-a],y2);
sum(y)
%std(y)
x_     = circshift(x,+1);
% y = r*[1,-r,1] y + (1-r^2)x
% ((1-2r+r^2) I + r [1,-2,1]))y + (1-r^2)x
% (I + r/(1-2r+r^2) D) y + (1-r^2)/(1-2*r+r^2)x

y(:,3) = filter(b*b,[-a,1+a^2,-a],x);
% note that the shifted two-sided filter is not stable, but the inverse filter is,
% so we test it by inversion
y_ = circshift(y(:,2),-1);
x_ = filter([-a,1+a^2,-a],(1-2*a+a^2),(1-2*a+a^2)/(1-a^2)*y_);
x_ = filter([0,1-2*a+a^2,0]-a*[1,-2,1],(1-2*a+a^2),(1-2*a+a^2)/(1-a^2)*y_);
a_ = a/(1-2*a+a^2);
x_ = filter([0,1,0]-a_*[1,-2,1],1,y_);
sum(x_)
pause
%y(:,3) = circshift(y(:,2),0);
%y = y./std(y);
%y(:,2) = filter([0,1/a],[1,-(1+a^2)/a,1],x); %,[0.2,0.1]);
%y(:,3) = filter([1],[1,-2*a,a^2],x);
%y = y./std(y);

subplot(2,3,4)
%parcorr(y);
%plot(y)
plot([x,x_])
%axis(ax);

end

if (0)

% two sided highass
yhp  = highpass1d_implicit(x,a);
y    = (x-filter(sqrt(1-a^2),[1,-a],x));
%yhp_ = filter([1-sqrt(1-a^2),-a]b*b*[a, -a^2, a],[-a, 1+a^2,-a],x);

subplot(2,3,5)
parcorr(yhp);

% two sided lowpass
y=bandpass1d_implicit(x,a);
subplot(2,3,6)
parcorr(y);
end

end
