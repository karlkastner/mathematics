a = 2;
b = 3;

fx = linspace(0,10)';
S = gampdf(fx,a,b);
[fc,Sc] = gamma_mode(a,b);

clf
plot(fx,S);
hold on
plot(fc,Sc,'*')


