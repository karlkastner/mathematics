% Thu 16 Feb 10:21:21 CET 2017
% Karl Kastner, Berlin
win_C = {'rectwin','triwin','hanwin','lanczoswin'}

for idx=1:5
	figure(idx)
	clf
end

n=2^10;
x = 10*(0:n-1)'/n;
[f T mask] = fourier_axis(x);
p=linspace(0,1,3);
y = [];
for idx=1:length(p);
	%T0=p(idx)*T(10)+(1-p(idx))*T(11);
	%T0 = exp( p(idx)*log(T(10)) + (1-p(idx))*log(T(11)));
	T0 = wharmean([p(idx); (1-p(idx))], T(10:11))
	y(:,idx) = sin(2*pi*x/T0);
end

dr =[];
r =[];
for idx=1:length(win_C)
win = n*feval(win_C{idx},x);
wy = bsxfun(@times,win,y);
F = abs(fft(wy));

figure(1);
subplot(2,2,idx)
semilogy(T(mask),F(mask,:));
title(win_C{idx});

figure(3);
subplot(2,2,1)
semilogy(T(mask),2/n*F(mask,2));
hold on

figure(2);
subplot(2,2,idx)
plot(2*F(mask,:)/n)
xlim([5 15])
ylim([0 1]);
title(win_C{idx});

figure(3);
subplot(2,2,2)
plot(2*F(mask,1)/n);
xlim([5 15])
hold on
legend(win_C{:});

dr(idx) = max(F(:,2))/min(F(:,2));
[maxF mdx] = max(F(:,1));
r(idx)  = maxF/F(mdx+1,1);

end

figure(4);
clf();
subplot(2,2,1);
bar(dr)
set(gca,'yscale','log');
xticklabel(win_C{:});
title('Dynamic range');

subplot(2,2,2);
bar(r);
set(gca,'yscale','log')
xticklabel(win_C{:});
title('Resolution');

%
% increased separation with increased sampling frequency
%
N = 2.^(5:7);
p = 0.5;
y = [];
T0 = wharmean([0.5 0.5],[10/9,1])
figure(5);
clf
dr = [];
for idx=1:length(N)
	n = N(idx);
	x = 10*(0:n-1)'/n;
	[f T mask] = fourier_axis(x);
	y = sin(2*pi*x/T0);
	F = 2/n*abs(fft(y));
	subplot(2,2,1)
	semilogy(T,F);
	hold on
	dr(idx) = max(F)/min(F);
end
subplot(2,2,2)
loglog(N,dr,'o')

%
% increased separation with increased length
%
%N = 2.^(5:7);
L = 2.^(0:0.1:18); %10*2.^(-2:5)';
y = [];
p = 0.5;
T0 = wharmean([p,1-p],[10,11])/10
%T0 = wharmean([0.5 0.5],[10/9,1])
figure(6);
clf
dr = [];
n0 = 4;
for idx=1:length(L)
	n = n0*L(idx)/L(1);	
	x = L(idx)*(0:n-1)'/n;
	[f T mask] = fourier_axis(x);
	y = sin(2*pi*x/T0);
	F = 2/n*abs(fft(y));
	if (1==mod(idx,10))
	subplot(2,2,1);
	semilogy(T,F);
	hold on;
	end
	[mF mdx] = max(F);
	mdx
	T(mdx)
	dr(idx) = max(F)/min(F);
end
xlim([0 10])
subplot(2,2,2)
loglog(L,dr,'.-')

figure(7);
clf
x=linspace(0,100)';
w = [];
for idx=1:length(win_C)
	w(:,idx) = feval(win_C{idx},x); 	
end
semilogy(abs(fftshift(fft(w,1024))))

