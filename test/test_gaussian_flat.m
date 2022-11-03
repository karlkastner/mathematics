 n = 1e6; L = 10; s=1; x = linspace(0,L,n)'; y = s*randn(n,1); S = periodogram(y,L,n); mean(S), std(S), Se=2*L/n, y = innerspace(0,20,1e5)'; 1/4*s*chi2pdf(2,y)'*y*(y(2)-y(1))*Se
>> fx=fourier_axis(x); df=fx(2)-fx(1); sum(S(fx>0))*df %df=1/L; n/2*Se*df = 1; 

ans =

    1.0000

>> fx=fourier_axis(x); df=fx(2)-fx(1); sum(S(fx>0))*df, df=1/L; Se = % n/2*Se*df = 1; 

