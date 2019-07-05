% 2017-07-07 11:21:26.936564308 +0200

% x = linspace(0,10,1e3)';
 x = linspace(0,10,1e3)';
 xi=linspace(-2,12,1e3);

 y = sin(3*x);
 y(x>4 & x < 6) = NaN;
 y(x>7 & x < 7.5) = NaN;
 y(x<1.1) = NaN;
 y(randi(length(y),100,1)) = NaN;
% y(:,2) = fixnan(x,y,1/10);
% y(:,2) = fixnan2(x,y,1/10);
 yi = interp1_limited(x,y,xi,1);
clf
 plot(xi,yi,'.');
hold on
 plot(x,y,'.');


