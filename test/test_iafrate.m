% Mon  4 Jul 15:14:42 CEST 2022
%
%
% note: the probability of being zero at xmin or having a zero crossing is then
% P(y(xmin)<0)*(1-Pcross) + Pcross

dx = 0.001;
L  = 10;
mu = 1;
m  = 1e4;
sx = 1;

x = (0:dx:L)';
n = length(x);
e = randn(n,m);
B = mu*x + sx*sqrt(dx)*cumsum(e);

xmin = 2;
B_ = B;
mdx=round(xmin/dx);
%[fdx,mdx] = min((x-xmin).^2);
%B_(x<xmin,:) = 1;
B_ = B_ - 0*B_(mdx,:);
B_(1:mdx-1,:)=1;
P = mean(cummin(B_)<0,2);
% & x>xmin,2);
% detect sign change
P_ = inner2outer(B(1:end-1,:).*B(2:end,:) < 0);
P_(1:mdx-1,:) = 0;
P_ = cummax(P_);
P_ = mean(P_,2);

P(:,2) = iafrate(xmin,x,mu,sx);
%P(:,3) = sqrt(P(:,2));
close all
plot(x,[P,P_])

