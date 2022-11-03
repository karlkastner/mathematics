 f=[10:10:40]; n=1e3; x = linspace(0,1,n)'; y = sum(cos(2*pi*x*f),2); [ff,FF] = bartlett_angle(y,1,10); plot(angle(FF(:,1:3)))
