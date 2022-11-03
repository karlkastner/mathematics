ne = 100;
n = 1e3;
x = randn(n,1)
y = randn(n,1)
xf = x;
yf = y;

clf
mu = mean([xf,yf]);                                     
                C  = cov([xf,yf]);                                      
                c  = [C(1),C(2),C(4),mu(1),mu(2)];                              
                [X, Y] = ellipse(c,[],ne);                                      
                %plot(10.^X,10.^Y,'--','linewidth',2,'color',cm(idx,:)); 
                plot(X,Y,'--','linewidth',2,'color',cm(idx,:)); 
axis equal
c=[1,2,3];x=randn(1e7,2); mean(hypot(x(:,1),x(:,2))<sqrt(-2*log((2*normcdf(-c))))), normcdf(c)-normcdf(-c)
n=1e5; x = randn(n,1); x(:,2) = randn(n,1) + 10*x; clf; plot(x(:,1),x(:,2),'.'), C = 1.52^2*cov(x); [v,e] = eig(C); e=e; x = x*v*inv(sqrt(e)); plot(x(:,1),x(:,2),'.'); axis equal; quantile(x,[0.16,0.84]), mean(hypot(x(:,1),x(:,2))<1)

