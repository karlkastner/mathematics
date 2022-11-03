% 2017-02-09 12:17:44.029177108 +0100

if(0)
n=100;
x = (-n:n)';
x=x/(2*n+1);

win = lanczoswin(x)
w2 = hanwin(x,0,0.2197);
subplot(2,2,1)
cla
plot(win);
hold on
plot(w2)
sum(win)
%[W(1,:)' A])

n = 2.^(ceil(log2(length(win)))+4);
subplot(2,2,2);
cla
f1 = abs(fft(win,n));
semilogx(f1);
hold on
f2=abs(fft(w2,n));
semilogx(f2);

z1 = fzero( @(x) interp1(1:n,f1,x,'linear')-1/sqrt(2), 30)
z2 = fzero( @(x) interp1(1:n,f2,x,'linear')-1/sqrt(2), 30)
z1/z2

% TODO-> equivalent filter length can be found from dof (kish)
%	-> but how can this be related to the cutoff frequency?
%	-> smoother windows have less sharp transition, so fc ~ 1/n_effective
 n=100; m=1e6; x = randn(n,m); w=lanczoswin(1:n); xl =(w*x); xm=mean(x); xh=hanwin(1:n)*x; 1./var(xl), 1./var(xh), 1./var(xm)
end

w_C = {
'rectwin'
'triwin'
'hanwin'
'kaiserwin'
'lanczoswin'
}
n = 1e6;
for idx=1:size(w_C,1)
	w_C{idx,2} = wdof(feval(w_C{idx},1:n))/n
end

