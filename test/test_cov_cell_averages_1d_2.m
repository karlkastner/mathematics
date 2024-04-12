% Tue 28 Nov 15:33:48 CET 2023
% first test, simple variance (x=0,y=0)
r = 0.35;
cfun = @(x,y) exp(-r*abs(x));
dx = 1;
m = 20;
x0 = 0:10;

ca = cov_cell_averages_1d(cfun,x0,dx,m,true)'
%cov_cell_averages_1d(cfun,x0,dx,m,false)

% manually
x = innerspace(-dx/2,dx/2,m);
x = flat(x);
for idx=1:length(x0)
C = cfun(x-x'+x0(idx));
c(idx) = mean(C,'all');
end
c


%if (0)
mu = 0;
sd = 1;
theta = 1./r;
nx = 100;
L = 100;
dx = 1;
[val,c] = ar1_1d_grid_cell_averaged_generate(mu,sd,theta,L,nx,m);
sd
c(1:11)'
c = c(1:11)'
c./ca
%end


