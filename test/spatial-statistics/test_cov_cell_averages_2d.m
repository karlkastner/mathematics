% Wed 29 Mar 21:46:40 CEST 2023
% Karl Kastner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% uni-directional
order = 9;
Lr = 1;
cfun = @(x,y) exp(-abs(x)/Lr)
%cfun = @(x,y) 0.5*(abs(x)<0.5);
L_ = 4;
%cfun = @(x,y) (1-abs(x/L_)).*(abs(x/L_)<1)
%cfun = @(x,y) 0;
dx = 1;
dy = dx;


x = [0,1,2,3,4,5]*dx;
y = [0,0,0,4,5];
c = cov_cell_averages_2d(cfun,x,y,dx,dy,order)
r=c/c(1)

n  = 1e5;
m  = 10;
x_ = dx*5*fourier_axis(0.5*m^2,m^2+1)';

%x_ = flat(innerspace(0,1,m)' + (-10:10)*Lr)';
m_ = length(x_);
c_ = cfun(x_,0);
% random
z = 1.*randn(n,m_);
%z(:,1) = 1;
T = sqrt(fft(c_,[],2));
e = ifft(T.*fft(z,[],2),[],2);
e = mid(e')';
m_ = m_ - 1;

c0 = corr(e);
% average over grid cells
%e =
ebar = zeros(n,m_/m);
xbar = zeros(1,m_/m);
xbar_ = zeros(1,m_/m);
cbar = zeros(1,m_/m);
x__ = fftshift(x_);
c__ = fftshift(c_);
for idx=1:(m_/m)
	for jdx=1:m
		xbar(idx) = xbar(idx) + x__((idx-1)*m+jdx);
		cbar(idx) = cbar(idx) + c_((idx-1)*m+jdx);
		xbar_(idx) = xbar_(idx) + x_((idx-1)*m+jdx);
		ebar(:,idx) = ebar(:,idx) + e(:,(idx-1)*m+jdx);
	end
end
xc=(0:(m_/m)-1)/(m_/m)*range(x_);
cbar = cbar/m;
ebar = ebar/m;
xbar = xbar/m;
xbar_ = xbar_/m;
c__ = cov(ebar)
r__ = c__./diag(c__);

subplot(2,2,1);
cla
plot(mid(x_),c0(1,:));
hold on
plot(mid(x_),c0_(1,:));
hold on
plot(x_,c_);
title('unaveraged');

subplot(2,2,2)
cla
plot(xc,c__(1,:),'o')
hold on
plot(x,c,'*')
xlabel('lag distance')
ylabel('cov');
legend('analytic','mc');
title('averaged');

subplot(2,2,3)
cla
plot(xc,r__(1,:),'o')
hold on
plot(x,r,'*')
ylabel('corr')
xlabel('lag distance')

