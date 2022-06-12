% 2022-04-09 14:42:42.274759731 +0200
I=[];
syms x p L;
L = 1;
for idx=1:6;
	I(idx,1) = int((1+cos(2*pi*x/L))^idx,x,0,L);
	I(idx,2) = nchoosek(2*idx,idx)/2^idx;
	I(idx,3) = (factorial(2*idx)/(factorial(idx)^2))/2^idx;
	I(idx,4) = (gamma(2*idx+1)./(gamma(idx+1).^2))./2.^idx;
	I(idx,5) = moments_fourier_power(idx,L);
end
I
clf
subplot(2,2,1)
plot(1:size(I,1),I,'*');
hold on
x = innerspace(sqrt(eps),7);
I = (gamma(2*x+1)./(gamma(x+1).^2))./2.^x;
plot(x,I)

% E(a*x)   = a*mu
% Var(a*x) = E(a*x - a*mu)^2 = a^2*x^2 - a^2*mu^2
p = linspace(0,10)';
mu = mean_fourier_power(p);
sd = std_fourier_power(p);
pause
subplot(2,2,2)
plot(p,[mu,sd])
subplot(2,2,3)
plot(p,[sd./mu])
%subplot(2,2,4)
%plot(p,[sd./mu].^2)
hold on
x = innerspace(0,1);
y = (1 + cos(2*pi*x));
for p=0:6
	s = std(y.^p,1);
	m = mean(y.^p);
	plot(p,s/m,'*');
	%, std_fourier_power(p)
end

%fac(4*idx)/fac(2*idx)^2 / (fac(2*idx)^2/fac(idx)^2) - 1
%fac(4*idx)*fac(idx)^2/fac(2*idx)^2
