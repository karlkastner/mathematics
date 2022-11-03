 n=100; x = innerspace(0,1)'; y = sin(10*pi*x); clf; plot(acf(y,99)); s = std(y); m=mean(y); y = y-mean(y); for k=1:n; a(k) = 1/(s^2*(n-k+1))*sum(y(1:n-(k-1)).*y(1+(k-1):n)); end; hold on; plot(a)

