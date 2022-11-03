n=100;
e = 0;
t = linspace(0,2*pi,n)';
s = sin(t)+e*randn(n,1);
c = cos(t)+e*randn(n,1);
sl = [];
R2 = [];
yp = [];
for idx=1:length(t)
	[sl(idx,:),intercept(idx,:),R2(idx,:),yp(idx,:)] = ls_perpendicular_offset([c(idx),0],[s(idx),0]);
end
figure(1)
subplot(2,2,1)
plot(t,[atan(sl),atan(c./s)])
ylabel('atan(slope)');

subplot(2,2,2)
plot(t,intercept)
ylabel('intercept');

subplot(2,2,3)
plot(t,R2)
ylabel('rmse2');

subplot(2,2,4)
plot(c,yp)
ylim([-1,1]); axis equal
%plot)

figure(2)
clf
t_ = atan(sl);
for idx=1:10:100
	subplot(3,4,(idx-1)/10+1)
	%plot(c(idx),s(idx),'*');
	plot([0,c(idx)],[0,s(idx)],'-*');
	hold on
	set(gca,'colororderindex',idx)
	%plot(cos(t_(idx)),sin(t_(idx)),'o')
	plot(cos(t_(idx)),sin(t_(idx)),'o')
	plot(-cos(t_(idx)),-sin(t_(idx)),'o')
	%plot(cos(sl(idx)),sin(sl(idx)),'o')
	ylim([-1.1,1.1]);
	axis equal
end
