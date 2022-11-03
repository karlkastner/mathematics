% Mon 25 Apr 13:56:34 CEST 2022
%
% requirement : impulse response = rho^|k| (autocorrelation)
% filtering: y = A^-1*x
% so the inverse A*y = (I - r*D2)y = x = [0 ... 1 ... 0] (impulse)
%
% therfore (D2/y).y = const = y, except at the origin

%  syms x y L positive; f = exp(-sqrt(x^2 + y^2)/L), r=simplify((diff(f,x,2)+diff(f,y,2))./f,'ignoreanalyticconstraints',true) 

%
% in 2d, this is not the case any more, even with exact derivatives:
% z = exp(-sqrt(x^2 + y^2)
% (D2/z).z = 1/L^2 - 1/(radius*L)
%
% i=1e1; j=1e1; r = 0.9; d2 = (r.^hypot(i-1,j) + r.^hypot(i+1,j) + r.^hypot(i,j-1) + r.^hypot(i,j+1) - 4*r.^hypot(i,j))./r.^(hypot(i,j)), (r.^(i-1) + r.^(i+1) - 2*r.^i)./r.^i  
%
% derivative must be defined along radius !!! -> which is strange, as the filter is unidirectional but the direction differs for each point
% d^2f/dR^2 + 1/R df/dR + 1/R^2 d2f/dtheta^2
% when derivative along theta not scaled : has different effect for different distance, nope, bc df/dtheta is zero
% when derivative along theta is left out : no smoothing between neigbouring rays
% -> maybe only the mixed derivative has to be left out

rho = 0.8;
n = 200;
x = zeros(2*n+1,1);
x(n+1,1) = 1;
y=lowpass1d_implicit(x,rho,1,true);
sum(y)
y=y/max(y(:));
r=rho./(1 - 2*rho + rho.*rho), k=(0:n)';
clf;

subplot(2,2,1)
plot(k,y(n+1:end,1));
hold on;
plot(k,rho.^k,'.');
xlim([0,20]);

D2 = derivative_matrix_2_1d(2*n+1,2*n,[],'circular');
subplot(2,3,4)
% (I - r*D)*y = 0
plot((D2*y)./y)



x = zeros(2*n+1,2*n+1);
x(n+1,n+1) = 1;
mode = 'fourier';
y=lowpass2d_implicit(x,rho,[],1,mode);
sum(y(:))
y=y/max(y(:));


subplot(2,3,6)
semilogy(k,y(n+1:end,n+1));
rho_ = y(n+1,n+2);
hold on;
plot(k,rho.^k,'.');
plot(k,rho_.^k,'.');
xlim([0,40]) 

subplot(2,2,4)
s=[]; R = linspace(0.011,0.999); k=-n:n; h=hypot(k,k'); for idx=1:length(R); s(idx,1) = sum(sum(R(idx).^h)); s(idx,2) = 2*pi/log(R(idx))^2; end; plot(R,s./s(:,2));% hline(1);

