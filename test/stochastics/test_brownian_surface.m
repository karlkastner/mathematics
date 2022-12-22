% c.f. schmidt 2015
% test 12.4.3 
% This is wiener, not brownian!

% note : 7.6 in powel and sharlow (covariance, variance seems wrong)

% alpha =
% covariance (cf. schmidt) 
% 	cov = 1/2 ( |s|^a + |t|^a + |s-t|^a )
% 	    = 1/2 (|s|^2H + |t|^2H + |s-t|^2H )
% hurst parameter
%	H = a/2 

% schmidt 397
%Xcurve : generated with newshams methods with covariance
%	c0 + c2|s - t|^2 - |s-t|^alpha
% correct
%Xtilde = Xcurve  - Xcurve0 + sqrt(2*c2)*t'*Z
% sum in x and y-direction


n = 512;
m = 1e3;

n = 32;
m = 1e2;

%n = 1024;
%m = 1;

% for BM = 0.5
a = 0.01;

d2 = 0;
id = (1:n)';
id0 = [n/2,n/2];
xx = sqrt(2*n)*brownian_field(a,round(4*n),m);
size(xx)
for jdx=1:m
if (1)
% x = randn(n,n)/sqrt(0.5*n);
% x = cumsum(x,1);
% x = cumsum(x,2);
	%x = sqrt(2*n)*brownian_field(0.5,round(4*n));
	x = xx(1:n,1:n,jdx);
else
	%
%	x = cumsum(randn(n,1))*cumsum(randn(1,n));
%	x  = x/sqrt(0.5*n);
%	x = x/(n^2);
end
 d2 = d2 + (x-x(id0(1),id0(2))).^2;
end
dij = hypot(id - id0(1),id'-id0(2));
d1 =  2*dij;
%abs(id-id0(1))+abs(id-id0(2))';

d2 = d2/m;
figure(1)
clf
subplot(2,2,1)
imagesc(x)
axis equal
subplot(2,2,2)
imagesc(d2)
axis equal
subplot(2,2,3)
imagesc(d1)

%d1 = abs(id - id0(1))*abs(id - id0(2))';
%$d1 = min(id-id0(1),id0(1)).*min(id-id0(2),id0(2))';
% covariance = prod min(x1,x2)
% variance, x1=x2,  = x*y


figure(2)

clf
subplot(2,2,1)
%imagesc(d2)
% diagonally
plot(diag(d2));
hold on
plot(diag(rot90(d2)))
plot(diag(d1));
%plot( abs(id-id0(1)).^2)
%plot(sqrt(2)*abs(id-id0(2)))

subplot(2,2,2)
plot(d2(id0(1),:))
hold on
%plot(abs(id-id0(1)))
plot(d1(id0(1),:))

subplot(2,2,3)
plot(d2(:,id0(2)))
hold on
plot(d1(:,id0(2)));
%abs(id-id0(2)))

%if(0)
subplot(2,2,4)
d=hypot(id-id0(1),id'-id0(2))+1;
%d = abs(id-id0(1)).*abs(id'-id0(2));
v=accumarray(round(d(:)),d2(:),[],@mean);
v_=accumarray(round(d(:)),d1(:),[],@mean);
plot(0:length(v)-1,v,'.')
hold on
plot(0:length(v)-1,v_,'.')
%plot((0:round(0.5*n))*sqrt(2))
%end

if (0)
clf
subplot(2,2,1)
d = (x-x(id0(1),:)).^2;
v = mean(d,2);
plot(v)
hold on
plot(abs(id-id0(1)));
subplot(2,2,2)
d = (x-x(:,id0(2))).^2;
v = mean(d,1);
plot(v)
hold on
plot(abs(id-id0(2)));

if (0)
subplot(2,2,3)
 d=hypot(id-id0(1),id'-id0(2))+1;
 v=1/n*accumarray(round(d(:)),(x(:)-x(id0(1),id0(2))).^2,[],@mean);
 plot([v,(1:length(v))'])
end
end
