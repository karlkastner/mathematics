% Tue  5 Dec 11:43:09 CET 2023
k=1.5;
 L=1000;
n=L^2;
 f = fourier_axis(L,n);
 [a,b] = logn_mode2par(1/k,1);
 S = lognpdf(f,a,b);
 [a,b]=logn_mode2par(k,1);
 P=lognpdf(f,a,b);
 plot(f,[S,P]);
df=1/L;
iS = sum(S)*df, 
iP = sum(P)*df, 
HD = 0.5*sum((sqrt(S)-sqrt(P)).^2)*df
HD=1-sum(sqrt(S.*P))*df 
HD = hellinger_distance(S,P,df)

