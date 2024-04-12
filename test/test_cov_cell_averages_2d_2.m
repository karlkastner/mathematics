% Tue 28 Nov 15:07:40 CET 2023
% first test, simple variance (x=0,y=0)
r = 0.2;
cfun = @(x,y) exp(-r*hypot(x,y));
dx = 1;
dy = 1;
m = 20;
n = m;
x0 = [0:3];
y0 = [0:3];
xx = repmat(cvec(x0),1,4);
yy = repmat(rvec(y0),4,1);
ca = cov_cell_averages_2d(cfun,xx,yy,dx,dy,m,n)

% manually
c = [];
for idx=1:length(x0)
for jdx=1:length(y0)
x = innerspace(-dx/2,dx/2,m);
x = repmat(cvec(x),1,n);
y = innerspace(-dy/2,dy/2,m);
y = repmat(rvec(y),m,1);
x = flat(x);
y = flat(y);
C = cfun(x-x'+x0(idx),y-y'+y0(jdx));
c(idx,jdx) = mean(C,'all');
end
end
c

if (0)
mu = 0;
sd = 1;
theta = 1./r;
nx = 1;
L = 1;
dx = 1;
ar1_2d_grid_cell_averaged_generate(mu,sd,theta,L*[1,1],nx*[1,1],m,n)
end

mu = 0;
sd = 1;
theta = 1./r;
nx = 100;
L = 100;
dx = 1;
[val,cc] = ar1_2d_grid_cell_averaged_generate(mu,sd,theta,L*[1,1],nx*[1,1],m,m);
cc=cc(1:4,1:4)'
%c = c(1:11)'
cc./ca

