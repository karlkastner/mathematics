% Fri 27 Aug 10:33:05 CEST 2021

f=200
n=10000

m=10
x = innerspace(sqrt(eps),1,n)'
x_ = innerspace(sqrt(eps),1/m,n/m)'
fx = fourier_axis(x);
fx_ = fourier_axis(x_);
a=cos(2*pi*f*x)+0*cos(4*pi*f*x)+1e-3*randn(n,1);
b=cos(2*pi*f*x+1/2*pi)
[c,c2]=coherence(a,b,m,false);
[c_,c2_]=coherence(a,b,m,true);
 
mean(c)
n_ = round(n/m);
c2__=mscohere(a,b,ones(n_,1),0,n_)

figure(1)
clf
subplot(2,2,1)
plot([a,b])

subplot(2,2,2)
plot(fx_,c2)
hold on
plot(fx(1:length(c2_)),c2_)
plot(fx_(1:length(c2__)),c2__)
subplot(2,3,4)
plot(fx_,c2)
subplot(2,3,5)
plot(c2_)
subplot(2,3,6)
plot(fx_(1:lenght(c2__)),c2__)

[mv,mdx] = max(c2);
fx_(mdx)
fx(mdx)
