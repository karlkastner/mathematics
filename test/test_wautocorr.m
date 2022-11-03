% Fri 11 Aug 14:04:46 CEST 2017

n = 100;
m = 1000;
rho = 0.975;
x = randar1(1,rho,n,m);
x2 = randar1(1,rho,n,m);
%x = sin(2*pi*(0:n-1)/(1.5e4));
x=x+flipud(x2);
%x = (1:n)'*ones(1,m);

k = n;
w = rand(n,m);

a = [];
%a = 0;
%for idx=1:size(x,2)
%	%a(:,1) = mean(autocorr_man5(x,m-1),2);
%	a = a+autocorr(x(:,idx),k-1)/size(x,2);
%end
s2 = var(x);
s2_ = s2/sum(s2);;
a(:,1) = autocorr_man4(x,k,[],[],false)*s2_';
a(:,2) = autocorr_man4(x,k,[],[],true)*s2_';
%w = (1:n)';
%w = (1:n).^2';
%w = (1:n)'-n/3;
s2 = wvar(w,x);
s2_ = s2/sum(s2);;
a(:,3) = wautocorr(w,x,k,[],[],false)*s2_'
a(:,4) = wautocorr(w,x,k,[],[],true)*s2_';
%a(:,2) = mean(wautocorr(w,x,m,[],true),2);

figure(2)
clf()
subplot(2,2,1)
plot(0:k-2,a(1:end-1,:),'.')
a(end-2:end-1,:)
legend('a ub','a b', 'wa ub','wa b')
hold on
plot(0:k-2,acfar1(rho,n,0:k-2,true));
plot(0:k-2,acfar1(rho,n,0:k-2,false));

% full
r = true;
a(:,1) = mean(autocorr_man4(x,k,[],r,false),2);
a(:,2) = mean(autocorr_man4(x,k,[],r,true),2);
%w = (1:n)';
w = rand(n,m);
%w = (1:n).^2';
%w = (1:n)'-n/3;
a(:,3) = mean(wautocorr(w,x,k,[],r,false),2);
a(:,4) = mean(wautocorr(w,x,k,[],r,true),2);

subplot(2,2,2)
legend('a ub','a b', 'wa ub','wa b')
plot(0:k-2,a(1:end-1,:),'.')
a(end-2:end-1,:)
legend('a ub','a b', 'wa ub','wa b')

